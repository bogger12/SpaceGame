using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CelestialBody : MonoBehaviour {

    public float mass;



    public bool rotate = true;
    public bool hasOrbit;
    public bool hasOrbitalLine;

    public float rotateSpeed;

    public Vector3 initVelocity;
    public GameObject bodyOfInfluence = null;

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
        lr = gameObject.AddComponent<LineRenderer>() as LineRenderer;
        lr.loop = true;
        lr.startWidth = lineWidth;
        lr.endWidth = lineWidth;
        lr.startColor = lineColor;
        lr.endColor = lineColor;
        lr.material = new Material(Shader.Find("Sprites/Default"));
        lr.sortingLayerID = SortingLayer.NameToID("Background");
    }


    // Start is called before the first frame update
    void Start()
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
    }

    // Update is called once per frame
    void Update()
    {
        if (hasOrbit) {
            orbit.CalculatePositionVelocityatTime(Time.time);
            //transform.position = orbit.getPosition();
            transform.position = GameSystem.VPixelSnap(orbit.getPosition());

            orbit.DrawOrbitalLine(lineRenderer, 50);
        }

        if (rotate) GameSystem.Rotate(transform, rotateSpeed * Time.deltaTime);

    }
}
