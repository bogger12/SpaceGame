using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RocketScript1 : MonoBehaviour {
    public GameObject bodyOfInfluence;
    public Rigidbody2D myRigidbody2D;


    public GameObject temp1;
    public GameObject temp2;

    const float gravconst = 6.6725985E-11f; // fundamental universal constant

    // Orbital Elements of the Object's orbit

    Vector2 e_; // Eccentricity Vector

    float e; // Eccentricity

    float a; // Semi-Major Axis (m)

    // Don't need inclination
    // Don't need longditude of ascending node

    float w; // Longtitude of periapsis (rad)

    float T; // Orbital Period

    float t0; // Time of periapis (initial)




    // Start is called before the first frame update
    void Start() {
        //myRigidbody2D.AddForce(new Vector2(2f, 15f), ForceMode2D.Impulse);
        myRigidbody2D.velocity = new Vector2(0f, 15f);

    }

    // Update is called once per frame
    private void Update() {
        if (Input.GetKey(KeyCode.W)) {
            DoThrust();
        }
        if (Input.GetKey(KeyCode.S)) {
            DoNThrust();
        }

        //float timesinceperiapsis = Time.time-t0;
        getElipse(transform.position);
        //GetPositionVelocityatTime(timesinceperiapsis);
    }

    void FixedUpdate() {
        myRigidbody2D.AddForce(GetGravityForce(), ForceMode2D.Force);
        //myRigidbody2D.velocity += GetGravityForce()* 0.02f;

    }

    void DoThrust() {
        Vector3 v = myRigidbody2D.velocity;
        myRigidbody2D.AddForce(v.normalized / 10);
    }

    void DoNThrust() {
        Vector3 v = myRigidbody2D.velocity;
        myRigidbody2D.AddForce(-v.normalized / 10);
    }

    Vector2 GetGravityForce() { // Part of this is written incorrectly causing F2 to move around
        Vector2 gravityvector = (bodyOfInfluence.transform.position - transform.position); // Vector of direction from this entity to its gravitational influence
        float radius = gravityvector.magnitude; // covnerting km to m
        Debug.Log(string.Format("gravityvector: {0}", gravityvector));
        Debug.Log(string.Format("radius: {0}", radius));

        float bodymass = 5.9722E12f;

        float gravacc = (gravconst * bodymass * myRigidbody2D.mass) / (radius * radius);

        return myRigidbody2D.mass * gravacc * gravityvector.normalized;

        //Debug.Log(string.Format(
        //    "Acc of Gravity (m/s-2): {0}\n" +
        //    "Radius (m): {1}",
        //    gravacc, radius));

    }

    //bool isInElipse(Vector2 r) {
    //    Vector2 v = myRigidbody2D.velocity;


    //    //float sgp = myRigidbody2D.mass * 1000 * gravconst; // Standard gravitational parameter;
    //    float sgp = 5.9722E12f * gravconst; // Standard gravitational parameter;

    //    float a = (sgp * r.magnitude) / (2 * sgp - r.magnitude * (v * v).magnitude); // Semi-major axis

    //    //Vector2 h = v * r;
    //    float h = r.x * v.y - r.y * v.x;

    //    Vector2 e = (r / r.magnitude) - (v * h) / (sgp); // Eccentricity vector

    //    Vector2 F2 = -2 * a * e;

    //    Debug.Log(string.Format(
    //        "a: {0}\n" +
    //        "e: {1}" +
    //        "F2: {2}",
    //        a, e.normalized, F2));
    //    Debug.Log((F2 - r).magnitude + r.magnitude);
    //    Debug.Log(2 * a);
    //    Debug.Log(e.magnitude);


    //    Debug.DrawLine(transform.position, transform.position + (Vector3)v/1000, Color.red);
    //    //Debug.DrawLine(new Vector3(0,5,0), new Vector3(0, 10, 0), Color.red);
    //    temp1.transform.position = F2;
    //    temp2.transform.position = (Vector2)bodyOfInfluence.transform.position + a*e.normalized;

    //    return (F2 - r).magnitude + r.magnitude == 2 * a;
    //}

    void getElipse(Vector2 r) {
        // r_ = position vector, v_ = velocity vector
        Vector2 v = myRigidbody2D.velocity;

        float sgp = 5.9722E12f * gravconst; // Standard gravitational parameter;
        Debug.Log(string.Format("sgp: {0}", sgp));


        Vector2 e = ((v.magnitude * v.magnitude - sgp / r.magnitude) * r - Vector2.Dot(r, v) * v) / sgp;

        float specificmechanicalenergy = (v.magnitude * v.magnitude) / 2 - (sgp / r.magnitude);

        float a = -(sgp / (2 * (specificmechanicalenergy)));


        Debug.Log(string.Format("specific mechanical energy: {0}", specificmechanicalenergy));
        Debug.Log(string.Format("eccentricity: {0}", e.magnitude));
        Debug.Log(string.Format("semi-major axis: {0}", a));

        Vector3 F2 = -2 * a * e;
        temp1.transform.position = F2;

        Debug.DrawLine(transform.position, transform.position + (Vector3)v, Color.red);


        float timesinceperiapsis = 0; // unit: seconds - t
        // h_ = angular momentum
        Vector3 h = Vector3.Cross(r, v);

        Vector3 nhat = Vector3.Cross(Vector3.forward, h);

        float orbitalPeriod = 2 * Mathf.PI * Mathf.Sqrt((a * a * a) / sgp);

        float meanmotion = (2 * Mathf.PI) / orbitalPeriod; // unit: radians - n

        float M = timesinceperiapsis * meanmotion;

        float E = M; // Eccentric Anomaly, unit: radians - E
        while (true) {
            var dE = (E - e.magnitude * Mathf.Sin(E) - M) / (1 - e.magnitude * Mathf.Cos(E));
            E -= dE;
            if (Mathf.Abs(dE) < 1e-6) break;
        }

        float P = a * (Mathf.Cos(E) - e.magnitude);
        float Q = a * Mathf.Sin(E) * Mathf.Sqrt(1 - Mathf.Pow(e.magnitude, 2));





        float angle = 60;

        Vector2 drawfromangle = Vector2.zero;

        Vector2 lastpoint = Vector2.zero;

        for (int i = 0; i < 50; i += 5) {
            float seperationdistance = ((h.magnitude * h.magnitude) / (myRigidbody2D.mass * sgp)) * (1 / (1 + e.magnitude * Mathf.Cos(i)));
            Debug.Log("seperationdistance: " + seperationdistance);

            Vector2 lineend = new Vector2(Mathf.Cos(i) * seperationdistance, Mathf.Sin(i) * seperationdistance);
            Debug.DrawLine(lastpoint, lineend, Color.blue);
            lastpoint = lineend;

        }



    }

    // Calculate Orbital Elements from Position and Velocity vectors

    void CalculateOrbitalElements(Vector2 r_, Vector2 v_) {

        float r = r_.magnitude;
        float v = v_.magnitude;


        float sgp = 5.9722E12f * gravconst; // Standard gravitational parameter;

        this.e_ = ((v * v - sgp / r) * r_ - Vector2.Dot(r_, v_) * v_) / sgp;
        this.e = e_.magnitude;

        float En = (v * v) / 2 - (sgp / r); // Specific mechanical energy

        this.a = -(sgp / (2 * (En)));

        this.T = 2 * Mathf.PI * Mathf.Sqrt((a * a * a) / sgp);


    }

    void GetPositionVelocityatTime(float timeSincePeriapsis) {

        float n = (2 * Mathf.PI) / T; // Mean motion (rad)
        float M = n * (timeSincePeriapsis); // Mean Anomaly (rad)

        float E = M; // Eccentric Anomaly (rad) - E
        while (true) {
            var dE = (E - e * Mathf.Sin(E) - M) / (1 - e * Mathf.Cos(E));
            E -= dE;
            if (Mathf.Abs(dE) < 1e-6) break;
        }

        float P = a * (Mathf.Cos(E) - e);
        float Q = a * Mathf.Sin(E) * Mathf.Sqrt(1 - Mathf.Pow(e, 2));

        var vP = -a * Mathf.Sin(E) * n / (1 - e * Mathf.Cos(E));
        var vQ = a * Mathf.Cos(E) * Mathf.Sqrt(1 - e * e) * n / (1 - e * Mathf.Cos(E));


        Vector2 r_ = new Vector2(P, Q);
        Vector2 v_ = new Vector2(vP, vQ);

        transform.position = r_;
        //myRigidbody2D.velocity = 

    }


}
