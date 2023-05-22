using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CelestialBody : MonoBehaviour {


    [Range(0, 10)]
    public float rotateSpeed;
    public bool rotate = true;


    public bool hasOrbit;

    public float mass;

    // If it has an orbit:
    public Orbit orbit;

    public Vector3 initVelocity;
    public GameObject bodyOfInfluence = null;


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
        }
    }

    // Update is called once per frame
    void Update()
    {
        if (hasOrbit) {
            orbit.CalculatePositionVelocityatTime(Time.time);
            transform.position = orbit.getPosition();
        }

        //orbit.DrawOrbitalLine(lineRenderer, orbitalLineNumberOfPoints);

        if (rotate) GameSystem.Rotate(transform, rotateSpeed * Time.deltaTime);

    }
}
