using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UIElements;

public class Orbit {

    const float gravconst = 6.6725985E-11f; // fundamental universal constant


    private Transform bodyOfInfluence;

    private Vector3 position;
    private Vector3 velocity;
    private float mainBodyMass;
    private float bodyOfInfluenceMass;

    public bool landed = false;

    // Standard gravitational parameter - Mu
    private float sgp;


    // == Orbital Elements of the Object's orbit == 

    // Eccentricity Vector
    private Vector3 e_; 

    // Eccentricity
    private float e;

    // Semi-Major Axis (m)
    private float a;

    // Argument of periapsis (rad)
    private float w;

    // Orbital Period
    private float T;

    // Epoch (Initial time at periapsis) - I use this to store the last time the orbital elements were changed
    private float t0;

    // Inclination
    private float i;

    // Longtitude of Right Ascending Node (RAAN)
    private float omega;

    // True Anomaly (rad)
    private float f;

    // == Extra variables ==

    private Vector3 pe_;
    private Vector3 ap_;
    private float pe;
    private float ap;




    // Initialises the state variables in the orbit, given BoI, and calculates orbital elements
    public Orbit(Vector3 position, Vector3 velocity, Transform bodyOfInfluence, float mainBodyMass, float bodyOfInfluenceMass) {
        this.mainBodyMass = mainBodyMass;
        this.bodyOfInfluenceMass = bodyOfInfluenceMass;
        this.bodyOfInfluence = bodyOfInfluence;

        this.sgp = bodyOfInfluenceMass * gravconst;

        CalculateOrbitalElementsFromPositionVelocity(position, velocity);
    }

    // Constructor that gives us the 6 Orbital Elements + BOI Information
    public Orbit(float eccentricity, float semiMajorAxis, float argumentOfPeriapsis,
        float inclination, float longtitudeOfAscendingNode, float trueAnomaly,
        Transform bodyOfInfluence, float mainBodyMass, float bodyOfInfluenceMass
    ) {
        this.mainBodyMass = mainBodyMass;
        this.bodyOfInfluenceMass = bodyOfInfluenceMass;
        this.bodyOfInfluence = bodyOfInfluence;

        this.sgp = bodyOfInfluenceMass * gravconst;

        this.e = eccentricity;
        this.a = semiMajorAxis;
        this.w = argumentOfPeriapsis;
        this.i = inclination;
        this.omega = longtitudeOfAscendingNode;
        this.f = trueAnomaly;

        this.e_ = RotateFromOrbitalPlaneTo3D(e, 0f);


        this.T = 2 * Mathf.PI * Mathf.Sqrt((a * a * a) / sgp);

        float n = (2 * Mathf.PI) / T; // Mean motion (rad)

        if (f != 0) this.t0 = Time.time - CalculateMeanAnomalyFrom(f) / n;
        else this.t0 = 0;


        CalculatePositionVelocityatTime(Time.time);

        //CalculateOrbitalElementsFromPositionVelocity(position, velocity);
    }

    public Vector3 GetPosition() {
        return bodyOfInfluence.position + position;
    }
    public Vector3 GetLocalPosition() {
        return position;
    }
    public Vector3 GetVelocity() {
        //return bodyOfInfluence.velocity + velocity;
        return velocity;
    }

    // Calculates the Orbital Elements of the Orbit given a position and velocity
    public void CalculateOrbitalElementsFromPositionVelocity(Vector3 r_, Vector3 v_) {

        float r = r_.magnitude;
        float v = v_.magnitude;

        // Standard gravitational parameter;
        float sgp = bodyOfInfluenceMass * gravconst;

        // h = angular momentum
        Vector3 h = Vector3.Cross(r_, v_);

        // e = Eccentricity
        this.e_ = ((v * v - sgp / r) * r_ - Vector3.Dot(r_, v_) * v_) / sgp;
        this.e = e_.magnitude;

        if (e >= 1) {
            Debug.LogError(string.Format("The eccentricity is too high. Haven't added calculations for this yet. e = {0}", e));
        }

        // a = semi-major axis
        this.a = 1 / (2 / r - (v * v) / sgp);

        // T = Orbital Period
        this.T = 2 * Mathf.PI * Mathf.Sqrt((a * a * a) / sgp);

        // i = Inclination
        this.i = Mathf.Acos(h.z / h.magnitude);

        // (temp) nodevector - Node Vector - points towards the ascending node
        //Vector3 nodevector = Vector3.Cross(Vector3.forward, h);
        Vector3 nodevector = Vector3.right; // Hack for 2D

        // omega = Longtitude of Ascending Node
        this.omega = Mathf.Acos(nodevector.x / nodevector.magnitude);
        if (nodevector.y < 0) this.omega = 2 * Mathf.PI - omega;
        if (i == 0) {
            omega = 0;
            nodevector = Vector3.right;
            // Hack as there is no nodevector if h.z == 0
        }

        this.w = Mathf.Acos(Vector3.Dot(nodevector, e_) / (nodevector.magnitude * e));
        if (e_.z < 0) this.w = 2 * Mathf.PI - w;
        if (i == 0) {
            this.w = Mathf.Atan2(e_.y, e_.x);
            if (Vector3.Cross(r_, v_).z < 0) this.w = 2 * Mathf.PI - w;
            // Hack as there is no nodevector if h.z == 0
        }

        this.f = CalculateArgumentofPeriapsis(r_, v_);


        // Find M at this point. Then use that to find t0
        // M at t0 = 0. Thus the difference in t is M/n

        float M = CalculateMeanAnomalyFrom(f);

        float n = (2 * Mathf.PI) / T; // Mean motion (rad)

        t0 = Time.time - (M / n); // Sets the epoch to the time of periapsis

        if (GameSystem.LOG_ELEMENTS) {
            Debug.Log(ToString());
        }


    }

    // Calculates the position and velocity of mainBody given a time
    public Vector3 CalculatePositionVelocityatTime(float inputTime, bool setvariables=true, bool evenmotion=false) {

        float timeSincePeriapsis = inputTime - t0;

        float n = (2 * Mathf.PI) / T; // Mean motion (rad)
        float M = n * (timeSincePeriapsis); // Mean Anomaly (rad)

        float E = M; // Eccentric Anomaly (rad) - E
        if (!evenmotion) for (int i=0;i<1000;i++) {
            var dE = (E - e * Mathf.Sin(E) - M) / (1 - e * Mathf.Cos(E));
            E -= dE;
            if (Mathf.Abs(dE) < 1e-6) break;
        }

        // Calculate X, Y and Z values
        float P = a * (Mathf.Cos(E) - e);
        float Q = a * Mathf.Sin(E) * Mathf.Sqrt(1 - Mathf.Pow(e, 2));

        // Calculate vX, vY and vZ values
        float vP = -a * Mathf.Sin(E) * n / (1 - e * Mathf.Cos(E));
        float vQ = a * Mathf.Cos(E) * Mathf.Sqrt(1 - e * e) * n / (1 - e * Mathf.Cos(E));

        Vector3 r_ = RotateFromOrbitalPlaneTo3D(P, Q);
        Vector3 v_ = RotateFromOrbitalPlaneTo3D(vP, vQ);


        if (setvariables) {
            this.f = CalculateArgumentofPeriapsis(r_, v_);
        }

        if (setvariables) {
            position = r_;
            velocity = v_;
            return Vector3.zero;
        }
        else return r_;

    }

    public static Vector3 SimulateGravity(Vector3 pos, Vector3 BOIpos, float BOImass) {
        float sgp = gravconst * BOImass;
        Vector3 difference = BOIpos - pos;
        float r = difference.magnitude;
        Vector3 gravacc = Time.deltaTime * (sgp / (r * r)) * difference.normalized;
        return gravacc;
    }

    // Adds the acceleration of a given force and recalculates orbital elements - also adds step in position
    public void AddForce(Vector3 force) {
        Vector3 thrustacc = Time.deltaTime * force / mainBodyMass;
        // a = (GM)/r^2

        // The Vector3.zero would usually be the body of influence pos
        Vector3 difference = Vector3.zero - position;
        float r = difference.magnitude;
        Vector3 gravacc = Time.deltaTime * (sgp / (r * r)) * difference.normalized;

        velocity += gravacc;
        velocity += thrustacc;
        position += velocity * Time.deltaTime; // Adding acceleration of forces
        if (!landed) CalculateOrbitalElementsFromPositionVelocity(position, velocity);
    }

    // Draws the line of the orbit with the given number of vertices
    public void DrawOrbitalLine(LineRenderer lineRenderer, int numberOfPoints, bool pixelSnap) {
        lineRenderer.positionCount = numberOfPoints;

        //float n = (2 * Mathf.PI) / T; // Mean motion (rad)

        //float shipE = Mathf.Acos((e + Mathf.Cos(f)) / (1 + e * Mathf.Cos(f)));
        //if (f > Mathf.PI) shipE = 2 * Mathf.PI - shipE;

        //float reverseETimeSincePeriapsis = shipE / n;


        Vector3 position;

        for (int i = 0; i < numberOfPoints; i++) {
            position = bodyOfInfluence.position + CalculatePositionVelocityatTime(((T / (float)numberOfPoints) * (float)i + t0), false, true);
            if (pixelSnap) position = GameSystem.VPixelSnap(position);
            lineRenderer.SetPosition(i, position);
        }
    }

    Vector3 RotateFromOrbitalPlaneTo3D(float P, float Q) {

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

        return new Vector3(x, y, z);
    }

    float CalculateMeanAnomalyFrom(float f) {
        float M = Mathf.Atan2(-Mathf.Sqrt(1 - e * e) * Mathf.Sin(f), -e - Mathf.Cos(f)) + Mathf.PI - e * ((Mathf.Sqrt(1 - e * e) * Mathf.Sin(f)) / (1 + e * Mathf.Cos(f)));
        return M;
    }

    float CalculateArgumentofPeriapsis(Vector3 r_, Vector3 v_) {
        // We use argument of latitude if eccentricity is 0 (circular orbit)

        if (e != 0f) {
            float fcompute = Vector3.Dot(e_, r_) / (e * r_.magnitude);
            if (fcompute < -1f || fcompute > 1f) {
                Debug.LogError(string.Format("True anomaly fcompute is fucked btw (< -1); fcompute = {0}", fcompute));
            }
            float tempf = Mathf.Acos(fcompute);
            if (Vector3.Dot(r_, v_) < 0) tempf = 2 * Mathf.PI - tempf;

            return tempf;
        } else { // Circular orbit
            Vector3 nodevector = Vector3.right; // Hack for 2D

            float u = Mathf.Acos((Vector3.Dot(nodevector, r_) / (nodevector.magnitude * r_.magnitude)));

            return u;
        }


    }

    public void CalculateExtraVariables() {
        pe = (1 - e) * a;
        ap = (1 + e) * a;
        pe_ = RotateFromOrbitalPlaneTo3D(pe, 0);
        ap_ = RotateFromOrbitalPlaneTo3D(-ap, 0);
    }

    public override string ToString() {
        string outtext = "";
        outtext += string.Format("Position, Velocity: {0}, {1}\n", GetPosition(), GetVelocity());
        outtext += string.Format("eccentricity: {0}\n", e);
        outtext += string.Format("eccentricity vector: {0}\n", e_);
        outtext += string.Format("Semi-Major Axis: {0}m\n", a);
        outtext += string.Format("Period (s): {0}s\n", T);
        outtext += string.Format("Inclination (rad): {0}\n", i);
        outtext += string.Format("Longtitude of Ascending Node (rad): {0}\n", omega);
        outtext += string.Format("Argument of Periapsis (rad): {0}\n", w);
        outtext += string.Format("True Anomaly (rad): {0}\n", f);

        return outtext;
    }
}
