using System.Collections;
using System.Collections.Generic;
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

    //public bool hasSOI;
    public CircleCollider2D SphereOfInfluence;

    public float radiusSOI;
    public float RadiusSOI {
        get {
            return radiusSOI;
        }
        set {
            if (SphereOfInfluence != null) SphereOfInfluence.radius = value;
            radiusSOI = value;
        }
    }

    public bool rotate;
    public bool hasOrbit;
    public bool hasOrbitalLine;

    public float rotateSpeed;

    public Vector3 initVelocity;
    public GameObject bodyOfInfluence;

    public float lineWidth = 0.1f;
    public Color lineColor = Color.white;

    private Orbit orbit;
    private LineRenderer lineRenderer;

    public bool customOrbit;

    // If it has a custom orbit:
    public float e;
    public float a;
    public float w;
    public float i;
    public float omega;
    public float f;

    void SetupLineRenderer(ref LineRenderer lr) {
#if UNITY_EDITOR 
        if (lr == null) {
            lr = gameObject.AddComponent<LineRenderer>() as LineRenderer;
        }
#else
        lr = gameObject.AddComponent<LineRenderer>() as LineRenderer;
#endif
        lr.loop = true;
        SetLineWidth(lr, lineWidth);
        lr.startColor = lineColor;
        lr.endColor = lineColor;
        lr.material = new Material(Shader.Find("Sprites/Default"));
        lr.sortingLayerID = SortingLayer.NameToID("Background");
    }

    void SetLineWidth(LineRenderer lr, float lineWidth) {
        lr.startWidth = lineWidth;
        lr.endWidth = lineWidth;
    }

    void SetupSOI(ref CircleCollider2D SOI) {
#if UNITY_EDITOR
        if (SOI == null) {
            SOI = gameObject.AddComponent<CircleCollider2D>() as CircleCollider2D;
        }
#else
        SOI = gameObject.AddComponent<CircleCollider2D>() as CircleCollider2D;
#endif
        SOI.radius = radiusSOI;

    }

    // Start is called before the first frame update
    public void Start()
    {
        if (hasOrbit && !customOrbit) {
            orbit = new Orbit(
                transform.position,
                initVelocity,
                bodyOfInfluence.transform,
                mass,
                bodyOfInfluence.GetComponent<CelestialBody>().mass
            );
        } else if (hasOrbit && customOrbit) {
            orbit = new Orbit(e, a, w, i, omega, f, bodyOfInfluence.transform, mass, bodyOfInfluence.GetComponent<CelestialBody>().mass);
        }
        if (hasOrbitalLine) SetupLineRenderer(ref lineRenderer);

        SetupSOI(ref SphereOfInfluence);
    }

    // Update is called once per frame
    void Update()
    {
        if (hasOrbit) {
            orbit.CalculatePositionVelocityatTime(Time.time);
            //transform.position = orbit.getPosition();
            //transform.position = GameSystem.VPixelSnap(orbit.GetPosition());
            transform.position = orbit.GetPosition();

            orbit.DrawOrbitalLine(lineRenderer, 50, true);
        }

        if (rotate) GameSystem.Rotate(transform, rotateSpeed * Time.deltaTime);

        SetLineWidth(lineRenderer, GameSystem.pixelUnit * GameSystem.screenScale/2);
    }
}
