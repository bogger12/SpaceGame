using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;

public class CelestialBody : MonoBehaviour {

    public float Radius {
        get { return radius; }
        set {
            GetComponent<SpriteRenderer>().size = new Vector2(value * 2, value * 2);
            GetComponent<CircleCollider2D>().radius = value;
            radius = value;
            //Debug.Log(GetComponent<SpriteRenderer>().drawMode);
        }
    }

    public CircleCollider2D bodyCollider;

    public Vector3 initVelocity;
    public CelestialBody bodyOfInfluence;

    public float lineWidth = 0.1f;
    public Color lineColor = Color.white;

    private Orbit orbit = null;
    private LineRenderer lineRenderer;

    public bool rotate = true;
    public bool hasOrbit = true;
    public bool customOrbit = true;
    public bool hasOrbitalLine = true;

    public bool hasCustomSOI = false;
    public float radiusSOI;

    public float mass;
    public float radius;
    public float density;
    public float rotatePeriod;
    // If it has a custom orbit:
    public float e;
    public float a;
    public float w;
    public float i;
    public float omega;
    public float f;


    // Start is called before the first frame update
    public void Start()
    {
        if (hasOrbit) {
            if (customOrbit) orbit = new Orbit(e, a, w, i, omega, f, bodyOfInfluence, mass);
            else orbit = new Orbit(transform.position, initVelocity, bodyOfInfluence, mass);
            if (!hasCustomSOI) radiusSOI = CalculateSphereOfInfluence();
        }
        if (hasOrbitalLine) SetupLineRenderer(ref lineRenderer);
        bodyCollider = GetComponent<CircleCollider2D>();
        bodyCollider.radius = radius;
        GetComponent<SpriteRenderer>().size = new Vector2(radius * 2, radius * 2);
    }

    // Update is called once per frame
    public void OrbitalUpdate()
    {
        if (hasOrbit) {
            orbit.CalculatePositionVelocityatTime(GameSystem.CurrentTime());
            //transform.position = GameSystem.VPixelSnap(orbit.GetPosition());
            transform.position = orbit.GetPosition();

            GameSystem.SetLineWidth(lineRenderer, GameSystem.pixelUnit * GameSystem.screenScale);
            orbit.DrawOrbitalLine(lineRenderer, 50, true);
            if (GameSystem.DEBUG) Debug.DrawCircle(transform.position, radiusSOI, 20, Color.red, 0);
        }

        if (rotate) GameSystem.Rotate(transform, 2f*Mathf.PI*(GameSystem.DeltaTime()/rotatePeriod) );

    }

    void SetupLineRenderer(ref LineRenderer lr) {
#if UNITY_EDITOR 
        if (lr == null) {
            lr = gameObject.AddComponent<LineRenderer>() as LineRenderer;
        }
#else
        lr = gameObject.AddComponent<LineRenderer>() as LineRenderer;
#endif
        lr.loop = true;
        GameSystem.SetLineWidth(lr, lineWidth);
        lr.startColor = lineColor;
        lr.endColor = lineColor;
        lr.material = new Material(Shader.Find("Sprites/Default"));
        lr.sortingLayerID = SortingLayer.NameToID("Background");
    }

    public Orbit GetOrbit() { return orbit; }
    public float GetSphereOfInfluence() { return radiusSOI; }
    public float CalculateSphereOfInfluence() {
        if (hasCustomSOI) return radiusSOI;
        else return 0.9431f * a * Mathf.Pow(mass / bodyOfInfluence.mass, (2f / 5f));
    }

    public int GetHeirarchyLevel() {
        if (bodyOfInfluence is null) return 0;
        else return bodyOfInfluence.GetHeirarchyLevel() + 1;
    }

    public CelestialBody GetRootBody() {
        if (bodyOfInfluence is null) return this;
        else return bodyOfInfluence.GetRootBody();
    }
}
