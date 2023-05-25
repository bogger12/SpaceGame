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


    void setupLineRenderer(ref LineRenderer lr) {
        lr = gameObject.AddComponent<LineRenderer>() as LineRenderer;
        lr.loop = true;
        lr.startWidth = lineWidth;
        lr.endWidth = lineWidth;
        lr.startColor = lineColor;
        lr.endColor = lineColor;
        lr.material = new Material(Shader.Find("Sprites/Default"));
        lr.sortingLayerID = 1;
    }


    // Start is called before the first frame update
    void Start()
    {
        if (hasOrbit && bodyOfInfluence != null) {
            orbit = new Orbit(
                transform.position,
                initVelocity,
                bodyOfInfluence.transform,
                mass,
                bodyOfInfluence.GetComponent<CelestialBody>().mass
            );

            if (hasOrbitalLine) setupLineRenderer(ref lineRenderer);
        }
    }

    // Update is called once per frame
    void Update()
    {
        if (hasOrbit) {
            orbit.CalculatePositionVelocityatTime(Time.time);
            transform.position = orbit.getPosition();

            orbit.DrawOrbitalLine(lineRenderer, 50);
        }

        if (rotate) GameSystem.Rotate(transform, rotateSpeed * Time.deltaTime);

    }
}
