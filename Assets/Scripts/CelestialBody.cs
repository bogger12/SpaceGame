using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;

public class CelestialBody : MonoBehaviour {

    public float mass;

    public float radius;
    public float Radius {
        get { return radius; }
        set {
            GetComponent<SpriteRenderer>().size = new Vector2(value * 2, value * 2);
            GetComponent<CircleCollider2D>().radius = value;
            radius = value;
            //Debug.Log(GetComponent<SpriteRenderer>().drawMode);
        }
    }
    public bool hasCustomSOI;
    public float radiusSOI;

    //public bool hasSOI;
    public CircleCollider2D bodyCollider;

    public bool rotate;
    public bool hasOrbit;
    public bool hasOrbitalLine;

    public float rotateSpeed;

    public Vector3 initVelocity;
    public CelestialBody bodyOfInfluence;

    public float lineWidth = 0.1f;
    public Color lineColor = Color.white;

    private Orbit orbit = null;
    private LineRenderer lineRenderer;

    public bool customOrbit;

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
            radiusSOI = orbit.CalculateSphereOfInfluence();
        }
        if (hasOrbitalLine) SetupLineRenderer(ref lineRenderer);
        bodyCollider = GetComponent<CircleCollider2D>();
    }

    // Update is called once per frame
    void Update()
    {
        if (hasOrbit) {
            orbit.CalculatePositionVelocityatTime(Time.time);
            //transform.position = GameSystem.VPixelSnap(orbit.GetPosition());
            transform.position = orbit.GetPosition();

            GameSystem.SetLineWidth(lineRenderer, GameSystem.pixelUnit * GameSystem.screenScale);
            orbit.DrawOrbitalLine(lineRenderer, 50, true);
        }

        if (rotate) GameSystem.Rotate(transform, rotateSpeed * Time.deltaTime);

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

    public int GetHeirarchyLevel() {
        if (bodyOfInfluence is null) return 0;
        else return bodyOfInfluence.GetHeirarchyLevel() + 1;
    }

    public CelestialBody GetRootBody() {
        if (bodyOfInfluence is null) return this;
        else return bodyOfInfluence.GetRootBody();
    }
}
