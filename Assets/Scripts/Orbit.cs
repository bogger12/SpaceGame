using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UIElements;
public class Orbit {

    public enum OrbitType { Circular, Elliptical, Parabolic, Hyperbolic }

    const float gravconst = 6.6725985E-11f; // fundamental universal constant


    private Transform bodyOfInfluence;

    private Vector3 position;
    private Vector3 velocity;
    private float mainBodyMass;
    private float bodyOfInfluenceMass;

    public bool landed = false;

    public OrbitType orbitType;

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

    // Epoch - Always set to time at periapsis
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
        orbitType = GetOrbitType(e);

        float n = CalculateMeanMotion(orbitType);

        // T = Orbital Period
        if (orbitType == OrbitType.Circular || orbitType == OrbitType.Elliptical) {
            this.T = (2 * Mathf.PI) / n;
        }
        else this.T = float.PositiveInfinity;

        if (f != 0) this.t0 = Time.time - CalculateMeanAnomaly(orbitType) / n;
        else this.t0 = 0;
        //Debug.Log(M);

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

        // e = Eccentricity
        this.e_ = ((v * v - sgp / r) * r_ - Vector3.Dot(r_, v_) * v_) / sgp;
        this.e = e_.magnitude;
        orbitType = GetOrbitType(e);

        // h = angular momentum
        Vector3 h = Vector3.Cross(r_, v_);

        // i = Inclination
        this.i = Mathf.Acos(h.z / h.magnitude);

        // (temp) nodevector - Node Vector - points towards the ascending node
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

        // a = semi-major axis
        this.a = 1 / (2 / r - (v * v) / sgp);

        this.f = CalculateTrueAnomaly(r_, v_, orbitType);

        // Find M at this point. Then use that to find t0
        // M at t0 = 0. Thus the difference in t is M/n

        float M = CalculateMeanAnomaly(orbitType);

        float n = CalculateMeanMotion(orbitType);

        // T = Orbital Period
        if (orbitType == OrbitType.Circular || orbitType == OrbitType.Elliptical) {
            this.T = (2 * Mathf.PI) / n;
        }
        else this.T = float.PositiveInfinity;

        t0 = Time.time - (M / n); // Sets the epoch to the time of periapsis


        if (GameSystem.LOG_ELEMENTS) {
            Debug.Log(ToString());
        }


    }

    // Calculates the position and velocity of mainBody given a time
    public Vector3 CalculatePositionVelocityatTime(float inputTime, bool setvariables=true, bool evenmotion=false) {

        float n = CalculateMeanMotion(orbitType);

        Vector2 P = Vector2.zero, V = Vector2.zero;

        if (orbitType == OrbitType.Circular || orbitType == OrbitType.Elliptical) {
            float E = CalculateEccentricAnomaly(inputTime, orbitType, evenmotion);
            // Calculate X, Y and Z values
            P.x = a * (Mathf.Cos(E) - e);
            P.y = a * Mathf.Sin(E) * Mathf.Sqrt(1 - e * e);

            // Calculate vX, vY and vZ values
            V.x = -a * Mathf.Sin(E) * n / (1 - e * Mathf.Cos(E));
            V.y = a * Mathf.Cos(E) * Mathf.Sqrt(1 - e * e) * n / (1 - e * Mathf.Cos(E));
        } else if (orbitType == OrbitType.Parabolic) {

            float f = CalculateTrueAnomalyFromTime(inputTime, orbitType);

            float h = Mathf.Sqrt(sgp * (a * (e * e - 1)));

            float r = ((h * h) / sgp) * (1 / (1 + Mathf.Cos(f)));

            float v = Mathf.Sqrt((2 * sgp) / r);

            P = GameSystem.PolarToCartesian(r, f);
            //V = PolarToCartesian(v, phi + (Mathf.PI / 2));

        } else if (orbitType == OrbitType.Hyperbolic) {
            float f = CalculateTrueAnomalyFromTime(inputTime, orbitType);
            float l = a * (e * e - 1);
            float r = -(l / (1 + e * Mathf.Cos(f))); // idk why tf but it has to be negative

            // Flight path angle
            float phi = Mathf.Atan((e * Mathf.Sin(f)) / (1 + e * Mathf.Cos(f)));

            float v = Mathf.Sqrt(sgp * ((2 / r) - (1 / a)));

            float vangle = f + (Mathf.PI / 2) - phi;

            P = GameSystem.PolarToCartesian(r, f);
            V = GameSystem.PolarToCartesian(v, vangle);
        }

        Vector3 r_ = RotateFromOrbitalPlaneTo3D(P.x, P.y);
        Vector3 v_ = RotateFromOrbitalPlaneTo3D(V.x, V.y);

        if (setvariables) {
            this.f = CalculateTrueAnomalyFromTime(inputTime, orbitType);
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
        

        //float n = (2 * Mathf.PI) / T; // Mean motion (rad)

        //float shipE = Mathf.Acos((e + Mathf.Cos(f)) / (1 + e * Mathf.Cos(f)));
        //if (f > Mathf.PI) shipE = 2 * Mathf.PI - shipE;

        //float reverseETimeSincePeriapsis = shipE / n;


        Vector3 position;

        if (orbitType==OrbitType.Circular||orbitType==OrbitType.Elliptical) {
            lineRenderer.positionCount = numberOfPoints;

            for (int i = 0; i < numberOfPoints; i++) {
                position = bodyOfInfluence.position + CalculatePositionVelocityatTime(((T / (float)numberOfPoints) * (float)i + t0), false, true);
                if (pixelSnap) position = GameSystem.VPixelSnap(position);
                lineRenderer.SetPosition(i, position);
            }
            lineRenderer.loop = true;
        }
        else if (orbitType == OrbitType.Parabolic) {

        } else if (orbitType == OrbitType.Hyperbolic) {
            float timerange = 50; // seconds
            lineRenderer.positionCount = numberOfPoints-1;
            for (int i = 1; i < numberOfPoints/2; i++) {
                float timeslice = (timerange / (float)numberOfPoints) * (timerange / ((float)i)) - 2 * timerange * timerange / (numberOfPoints * numberOfPoints);
                position = bodyOfInfluence.position + CalculatePositionVelocityatTime((t0+timeslice), false, true);
                if (pixelSnap) position = GameSystem.VPixelSnap(position);
                lineRenderer.SetPosition(i-1, position);
            }
            for (int i = 0; i <= numberOfPoints/2; i++) {
                float timeslice = (timerange / (float)numberOfPoints) * (timerange / ((float)i)) - 2*timerange*timerange /(numberOfPoints*numberOfPoints);
                if (i == 0) timeslice = 0;
                position = bodyOfInfluence.position + CalculatePositionVelocityatTime((t0 - timeslice), false, true);
                if (pixelSnap) position = GameSystem.VPixelSnap(position);
                lineRenderer.SetPosition(numberOfPoints-i-1, position);
            }
            lineRenderer.loop = false;
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

    float CalculateMeanMotion(OrbitType orbitType) {
        if (orbitType==OrbitType.Circular||orbitType==OrbitType.Elliptical) {
            float n = Mathf.Sqrt(sgp / (a * a * a));
            return n;
        } else if (orbitType==OrbitType.Parabolic) {
            //float n = Mathf.Sqrt(sgp / (a * a * a));
            //float n = Mathf.Sqrt(sgp / -(a * a * a));
            float rp = -a * (e - 1);
            float n = 2 * Mathf.Sqrt(sgp / rp*rp*rp);
            return n;
        }
        else if (orbitType==OrbitType.Hyperbolic) {
            float n = Mathf.Sqrt(sgp / -(a * a * a));
            return n;
        }
        return 0;
    }

    float CalculateMeanAnomaly(OrbitType orbitType) {
        // This assumes f is given

        if (orbitType==OrbitType.Circular||orbitType==OrbitType.Elliptical) {
            float M = Mathf.Atan2(-Mathf.Sqrt(1 - e * e) * Mathf.Sin(f), -e - Mathf.Cos(f)) + Mathf.PI - e * ((Mathf.Sqrt(1 - e * e) * Mathf.Sin(f)) / (1 + e * Mathf.Cos(f)));
            return M;
        } else if (orbitType == OrbitType.Parabolic) {
            float l = a * (e * e - 1);
            float q = l / 2;
            float D = CalculateEccentricAnomaly(0, orbitType);

            float M = q * D + (D * D * D) / 6f;
            Debug.Log(D);
            Debug.Log(M);
            return M;
        } else if (orbitType==OrbitType.Hyperbolic) {
            float first = Mathf.Tan(f / 2) / Mathf.Sqrt(((e + 1) / (e - 1)));
            float F = 2 * System.MathF.Atanh(first);

            float M = (e * System.MathF.Sinh(F)) - F;
            return M;
        }
        return 0;
    }

    float CalculateEccentricAnomaly(float currentTime, OrbitType orbitType, bool evenmotion=false) {
        // Parabolic assumes f is given

        float timeSincePeriapsis = currentTime - t0;

        if (orbitType==OrbitType.Circular) {
            float n = CalculateMeanMotion(orbitType);
            float M = n * (timeSincePeriapsis); // Mean Anomaly (rad)
            return M;
        }
        else if (orbitType == OrbitType.Elliptical) {
            // Using Newton's approximation method
            float n = CalculateMeanMotion(orbitType);
            float M = n * (timeSincePeriapsis); // Mean Anomaly (rad)

            float E = M; // Eccentric Anomaly (rad) - E
            if (!evenmotion) {
                if (e > 0 && e < 1) {
                    for (int i = 0; i < 1000; i++) {
                        var dE = (E - e * Mathf.Sin(E) - M) / (1 - e * Mathf.Cos(E));
                        E -= dE;
                        if (Mathf.Abs(dE) < 1e-6) break;
                    }
                }
            }
            return E;
        } else if (orbitType == OrbitType.Parabolic) {
            //float f = CalculateTrueAnomalyFromTime(currentTime, orbitType);
            float l = a * (e * e - 1);
            float D = Mathf.Sqrt(l)*Mathf.Tan(f / 2);
            return D;

        } else if (orbitType == OrbitType.Hyperbolic) {
            // Using Newton's approximation method
            float n = CalculateMeanMotion(orbitType);
            float M = n * (timeSincePeriapsis); // Mean Anomaly (rad)

            float F = M; // Eccentric Anomaly (rad) - E
            for (int i = 0; i < 1000; i++) {
                var dF = (e * System.MathF.Sinh(F) - F - M) / (e * System.MathF.Cosh(F)-1);
                F -= dF;
                if (Mathf.Abs(dF) < 1e-6) break;
            }
            return F;
        }

        return 0;

    }

    float CalculateTrueAnomaly(Vector3 r_, Vector3 v_, OrbitType orbitType) {
        // We use argument of latitude if eccentricity is 0 (circular orbit)

        if (orbitType == OrbitType.Circular) { // Circular Orbit
            Vector3 nodevector = Vector3.right; // Hack for 2D

            float u = Mathf.Acos((Vector3.Dot(nodevector, r_) / (nodevector.magnitude * r_.magnitude)));
            if (r_.z<0) u = 2 * Mathf.PI - u;

            return u;
        }else {
            float fcompute = Vector3.Dot(e_, r_) / (e * r_.magnitude);
            if (fcompute < -1 || fcompute > 1) {
                Debug.LogError(string.Format("True anomaly fcompute is fucked btw (< -1); fcompute = {0}", fcompute));
            }
            float tempf = Mathf.Acos(fcompute);
            if (Vector3.Dot(r_, v_) < 0) tempf = 2 * Mathf.PI - tempf;

            //float h = Mathf.Sqrt(sgp * (a * (e * e - 1)));
            Vector3 h_ = Vector3.Cross(r_, v_);
            float h = h_.magnitude;

            float q = Vector3.Dot(r_, v_);
            float fy = (h * q) / (r_.magnitude * sgp);
            float fx = (h * h) / (r_.magnitude * sgp) - 1;
            float tempf2 = Mathf.Atan2(fy, fx);

            return tempf2;
        }

    }

    float CalculateTrueAnomalyFromTime(float currentTime, OrbitType orbitType) {
        // We use argument of latitude if eccentricity is 0 (circular orbit)

        if (orbitType == OrbitType.Circular) { // Circular Orbit
            float E = CalculateEccentricAnomaly(currentTime, orbitType, false);
            return E;
        }
        else if (orbitType == OrbitType.Elliptical) { // Eliptical Orbit
            float E = CalculateEccentricAnomaly(currentTime, orbitType, false);
            return CalculateTrueAnomalyFromE(E, orbitType);
        }
        else if (orbitType == OrbitType.Parabolic) { // Parabolic Orbit

            //// Use Barker's equation to calcluate true anomaly 
            //float rp = -a * (e - 1);
            float h = Mathf.Sqrt(sgp * (a * (e * e - 1)));
            float rp = (h * h) / (2 * sgp);

            float A = (3f / 2) * Mathf.Sqrt((sgp) / (2 * rp * rp * rp)) * (currentTime - t0);

            float B = Mathf.Pow((A + Mathf.Sqrt(A * A + 1)), (1 / 3f));

            float tempf = 2 * Mathf.Atan(B - (1 / B));

            Debug.Log(rp);
            Debug.Log(A);
            Debug.Log(B);
            Debug.Log(tempf);

            return tempf;

        }
        else if (orbitType == OrbitType.Hyperbolic) { // Hyperbolic Orbit
            float F = CalculateEccentricAnomaly(currentTime, orbitType, false);
            return CalculateTrueAnomalyFromE(F, orbitType);
        }

        return 0;

    }

    float CalculateTrueAnomalyFromE(float E, OrbitType orbitType) {
        if (orbitType==OrbitType.Elliptical) {
            float B = e / (1 + Mathf.Sqrt(1 - e * e));
            float tempf = E + 2 * Mathf.Atan((B * Mathf.Sin(E)) / (1 - B * Mathf.Cos(E)));
            return tempf;
        } else if (orbitType == OrbitType.Hyperbolic) {
            float tempf = 2*Mathf.Atan(Mathf.Sqrt((e + 1) / (e - 1)) * System.MathF.Tanh(E / 2));
            return tempf;
        }
        return 0;
    }

    public void CalculateExtraVariables() {
        if (orbitType==OrbitType.Elliptical) {
            pe = (1 - e) * a;
            ap = (1 + e) * a;
            pe_ = RotateFromOrbitalPlaneTo3D(pe, 0);
            ap_ = RotateFromOrbitalPlaneTo3D(-ap, 0);
        } else if (orbitType==OrbitType.Hyperbolic) {
            pe = -a * (e - 1);
            pe_ = RotateFromOrbitalPlaneTo3D(pe, 0);
        }
    }

    public OrbitType GetOrbitType(float e) {
        if (e == 0) return OrbitType.Circular;
        else if (e > 0 && e < 1) return OrbitType.Elliptical;
        else if (e == 1) return OrbitType.Parabolic;
        else if (e > 1) return OrbitType.Hyperbolic;
        else return 0;
    }

    public override string ToString() {
        string outtext = "";
        outtext += string.Format("Position, Velocity: {0}, {1}\n", GetLocalPosition(), GetVelocity());
        outtext += string.Format("Eccentricity: {0}\n", e);
        outtext += string.Format("Eccentricity vector: {0}\n", e_);
        outtext += string.Format("Orbit Type: {0}\n", orbitType.ToString());
        outtext += string.Format("Epoch: {0}\n", t0);
        outtext += string.Format("Semi-Major Axis: {0}m\n", a);
        outtext += string.Format("Period (s): {0}s\n", T);
        outtext += string.Format("Inclination (rad): {0}\n", i);
        outtext += string.Format("Longtitude of Ascending Node (rad): {0}\n", omega);
        outtext += string.Format("Argument of Periapsis (rad): {0}\n", w);
        outtext += string.Format("True Anomaly (rad): {0}\n", f);

        return outtext;
    }
}
