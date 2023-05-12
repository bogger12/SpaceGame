using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RocketScript : MonoBehaviour {

    public UIController UIController;
    public LineRenderer lineRenderer;
    public GameObject bodyOfInfluence;

    const float gravconst = 6.6725985E-11f; // fundamental universal constant
    const float bodyOfInfluenceMass = 5.9722E12f;

    float sgp = bodyOfInfluenceMass * gravconst; // Standard gravitational parameter - Mu


    // Orbital Elements of the Object's orbit

    Vector3 e_; // Eccentricity Vector

    float e; // Eccentricity

    float a; // Semi-Major Axis (m)

    // Don't need inclination
    // Don't need longditude of ascending node

    float w; // Argument of periapsis (rad)

    float T; // Orbital Period

    float t0; // Epoch (Initial time at periapsis) - I use this to store the last time the orbital elements were changed

    float i; // Inclination

    float omega; // Longtitude of Right Ascending Node (RAAN)

    float f; // True Anomaly


    public Vector3 velocity = new Vector3(-1f, 15f);

    float mass = 1;



    // Start is called before the first frame update
    void Start() {
        CalculateOrbitalElementsFromPositionVelocity(transform.position, velocity);
        this.t0 = Time.time;
    }

    // Update is called once per frame
    private void Update() {
        if (Input.GetKey(KeyCode.W)) {
            addForce(velocity.normalized);
        }
        if (Input.GetKey(KeyCode.S)) {
            addForce(-velocity.normalized);
        }
        if (Input.GetKey(KeyCode.T)) {
            // Test update
            CalculateOrbitalElementsFromPositionVelocity(transform.position, velocity);
        }


        //CalculateOrbitalElementsFromPositionVelocity(transform.position, velocity);
        float timesinceperiapsis = Time.time - t0;

        DrawOrbitalLine();
        GetPositionVelocityatTime(timesinceperiapsis);
    }

    void addForce(Vector3 force) {
        Vector3 thrustacc = Time.deltaTime * force / mass;


        // F = (GMm)/r^2
        // a = (GM)/r^2
        Vector3 difference = bodyOfInfluence.transform.position - transform.position;
        float r = difference.magnitude;
        Vector3 gravacc = Time.deltaTime * (sgp / (r * r)) * difference.normalized;

        velocity += gravacc;
        velocity += thrustacc;
        transform.position += velocity * Time.deltaTime; // Adding acceleration of forces
        CalculateOrbitalElementsFromPositionVelocity(transform.position, velocity);
    }

    void CalculateOrbitalElementsFromPositionVelocity(Vector3 r_, Vector3 v_) {

        float r = r_.magnitude;
        float v = v_.magnitude;

        float sgp = bodyOfInfluenceMass * gravconst; // Standard gravitational parameter;

        // h = angular momentum
        Vector3 h = Vector3.Cross(r_, v_);

        this.e_ = ((v * v - sgp / r) * r_ - Vector3.Dot(r_, v_) * v_) / sgp;
        this.e = e_.magnitude;

        //float specificmechanicalenergy = (v * v) / 2 - (sgp / r); // Specific mechanical energy
        //this.a = -(sgp / (2 * (specificmechanicalenergy)));
        this.a = 1 / (2 / r - (v * v) / sgp);


        this.T = 2 * Mathf.PI * Mathf.Sqrt((a * a * a) / sgp);

        this.i = Mathf.Acos(h.z / h.magnitude);
        this.i = 0; // TODO: FIX ONCE 2D WORKS

        Vector3 nodevector = Vector3.Cross(Vector3.forward, h);

        this.omega = Mathf.Acos(nodevector.x / nodevector.magnitude);
        if (nodevector.y < 0) this.omega = 2 * Mathf.PI - omega;
        if (i == 0) {
            omega = 0;
            nodevector = Vector3.right;
        }
        //Debug.Log(string.Format("(i == 0) = {0}", (i == 0)));

        this.w = Mathf.Acos(Vector3.Dot(nodevector, e_) / (nodevector.magnitude * e));
        if (e_.z < 0) this.w = 2 * Mathf.PI - w;
        if (i == 0) {
            this.w = Mathf.Atan2(e_.y, e_.x);
            if (Vector3.Cross(r_, v_).z < 0) this.w = 2 * Mathf.PI - w;
        }

        float value = Vector3.Dot(e_, r_) / (e * r_.magnitude);
        if (value < -1f) {
            Debug.LogError(string.Format("True anomaly value is fucked btw (< -1); value = {0}", value));
            value = -1f;
        }
        else if (value > 1f) {
            Debug.LogError(string.Format("True anomaly value is fucked btw (> 1); value = {0}", value));
            value = 1f;
        }
        this.f = Mathf.Acos(value);
        if (Vector3.Dot(r_, v_) < 0) this.f = 2 * Mathf.PI - f;


        // Find M at this point. Then use that to find t0
        // M at t0 = 0. Thus the difference in t is M/n

        float M = Mathf.Atan2(-Mathf.Sqrt(1 - e * e) * Mathf.Sin(f), -e - Mathf.Cos(f)) + Mathf.PI - e * ((Mathf.Sqrt(1 - e * e) * Mathf.Sin(f)) / (1 + e * Mathf.Cos(f)));

        float n = (2 * Mathf.PI) / T; // Mean motion (rad)

        t0 = Time.time - (M / n); // TODO: idk if this will work

        Debug.Log(string.Format("Position, Velocity: {0}, {1}", r_, v_));

        Debug.Log(string.Format("nodevector: {0}", nodevector));

        Debug.Log(string.Format("eccentricity: {0}", e));
        Debug.Log(string.Format("Semi-Major Axis: {0}m", a));
        Debug.Log(string.Format("Period (s): {0}s", T));
        Debug.Log(string.Format("Inclination (rad): {0}", i));
        Debug.Log(string.Format("Longtitude of Ascending Node (rad): {0}", omega));
        Debug.Log(string.Format("Argument of Periapsis: {0}", w));
        Debug.Log(string.Format("True Anomaly: {0}", f));

        SetDebugText(r_, v_, h, nodevector, e_, a, T, i, omega, w, f);


    }

    Vector3 GetPositionVelocityatTime(float timeSincePeriapsis, bool setvariables=true, bool evenmotion=false) {

        float n = (2 * Mathf.PI) / T; // Mean motion (rad)
        float M = n * (timeSincePeriapsis); // Mean Anomaly (rad)

        float E = M; // Eccentric Anomaly (rad) - E
        if (!evenmotion) for (int i=0;i<1000;i++) {
            var dE = (E - e * Mathf.Sin(E) - M) / (1 - e * Mathf.Cos(E));
            E -= dE;
            if (Mathf.Abs(dE) < 1e-6) break;
        }


        //float sinE = Mathf.Sqrt(1 - e * e) * Mathf.Sin(f);
        //float cosE = e + Mathf.Cos(f);
        //float initE2 = Mathf.Atan2(sinE, cosE);

        //float aside = e + Mathf.Cos(f);
        //float hside = 1 + e * Mathf.Cos(f);
        //float myE = Mathf.Acos(aside / hside);
        //if (f > Mathf.PI) myE = 2 * Mathf.PI - myE;

        //Debug.Log(string.Format("My Calc initE: {0}", initE % (2 * Mathf.PI)));
        //Debug.Log(string.Format("My Calc initE2: {0}", initE2 % (2 * Mathf.PI)));

        //float myE = Mathf.Acos()

        //Debug.Log(string.Format("My Calc E: {0}", myE % (2 * Mathf.PI)));
        //Debug.Log(string.Format("Answer: {0}", E % (2 * Mathf.PI)));

        //E += initE; // When this is active, t0 must be set to the current time when recalculating orbital elements



        // Calculate X, Y and Z values
        float P = a * (Mathf.Cos(E) - e);
        float Q = a * Mathf.Sin(E) * Mathf.Sqrt(1 - Mathf.Pow(e, 2));
        // rotate by argument of periapsis
        float x = Mathf.Cos(w) * P - Mathf.Sin(w) * Q;
        float y = Mathf.Sin(w) * P + Mathf.Cos(w) * Q;
        // rotate by inclination
        float z = Mathf.Sin(i) * y;
        y = Mathf.Cos(i) * y;
        // rotate by longitude of ascending node
        float xtemp = x;
        x = Mathf.Cos(omega) * xtemp - Mathf.Sin(omega) * y;
        y = Mathf.Sin(omega) * xtemp + Mathf.Cos(omega) * y;

        // Calculate vX, vY and vZ values
        float vP = -a * Mathf.Sin(E) * n / (1 - e * Mathf.Cos(E));
        float vQ = a * Mathf.Cos(E) * Mathf.Sqrt(1 - e * e) * n / (1 - e * Mathf.Cos(E));
        // rotate by argument of periapsis
        float vx = Mathf.Cos(w) * vP - Mathf.Sin(w) * vQ;
        float vy = Mathf.Sin(w) * vP + Mathf.Cos(w) * vQ;
        // rotate by inclination
        float vz = Mathf.Sin(i) * vy;
        vy = Mathf.Cos(i) * vy;
        // rotate by longitude of ascending node
        float vxtemp = vx;
        vx = Mathf.Cos(omega) * vxtemp - Mathf.Sin(omega) * vy;
        vy = Mathf.Sin(omega) * vxtemp + Mathf.Cos(omega) * vy;


        Vector3 r_ = new Vector3(x, y, z);
        Vector3 v_ = new Vector3(vx, vy, vz);

        if (setvariables) {
            transform.position = r_;
            velocity = v_;
            return Vector3.zero;
        }
        else return r_;

    }

    void DrawOrbitalLine() {
        int numpoints = 100;
        lineRenderer.positionCount = numpoints;
        for (int i = 0; i < numpoints; i++) {
            lineRenderer.SetPosition(i, GetPositionVelocityatTime((T / (float)numpoints) * (float)i, false, true));
        }
    }

    void SetDebugText(Vector3 r_, Vector3 v_, Vector3 h, Vector3 nodevector, Vector3 e_, float a, float T, float i, float omega, float w, float f) {
        string outtext = "";
        outtext += string.Format("Position, Velocity: {0}, {1}\n", r_, v_);
        outtext += string.Format("nodevector (n): {0}\n", nodevector);
        outtext += string.Format("angular momentum (h): {0}\n", h);
        outtext += string.Format("eccentricity: {0}\n", e);
        outtext += string.Format("Semi-Major Axis: {0}m\n", a);
        outtext += string.Format("Period (s): {0}s\n", T);
        outtext += string.Format("Inclination (rad): {0}\n", i);
        outtext += string.Format("Longtitude of Ascending Node (rad): {0}\n", omega);
        outtext += string.Format("Argument of Periapsis (rad): {0}\n", w);
        outtext += string.Format("True Anomaly (rad): {0}\n", f);

        UIController.SetText(outtext);
    }

}
