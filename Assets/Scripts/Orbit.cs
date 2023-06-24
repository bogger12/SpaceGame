using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UIElements;

public class Orbit {

    public enum OrbitType { Circular, Elliptical, Parabolic, Hyperbolic }

    const float gravconst = 6.6725985E-11f; // fundamental universal constant


    private CelestialBody bodyOfInfluence;

    private Vector3 position;
    private Vector3 velocity;
    private float mainBodyMass;
    private float bodyOfInfluenceMass;

    public bool inStaticOrbit = true;

    public OrbitType orbitType;

    // Standard gravitational parameter - Mu
    private float sgp;

    // Specific Angular Momentum - h (conserved in closed system)
    private float h;

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
    public Orbit(Vector3 position, Vector3 velocity, CelestialBody bodyOfInfluence, float mainBodyMass) {
        this.mainBodyMass = mainBodyMass;
        this.bodyOfInfluenceMass = bodyOfInfluence.mass;
        this.bodyOfInfluence = bodyOfInfluence;

        this.sgp = bodyOfInfluenceMass * gravconst;

        CalculateOrbitalElementsFromPositionVelocity(position, velocity);
    }

    // Constructor that gives us the 6 Orbital Elements + BOI Information
    public Orbit(float eccentricity, float semiMajorAxis, float argumentOfPeriapsis,
        float inclination, float longtitudeOfAscendingNode, float trueAnomaly,
        CelestialBody bodyOfInfluence, float mainBodyMass
    ) {
        this.mainBodyMass = mainBodyMass;
        this.bodyOfInfluenceMass = bodyOfInfluence.mass;
        this.bodyOfInfluence = bodyOfInfluence;

        this.sgp = bodyOfInfluenceMass * gravconst;

        this.e = eccentricity;
        this.a = semiMajorAxis;
        this.w = argumentOfPeriapsis;
        this.i = inclination;
        this.omega = longtitudeOfAscendingNode;
        this.f = trueAnomaly;

        this.e_ = RotateFromOrbitalPlaneTo3D(e, 0f);
        orbitType = FindOrbitType(e);

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

    public Vector3 GetPosition() { return bodyOfInfluence.transform.position + position; }
    public Vector3 GetLocalPosition() { return position; }
    public Vector3 GetVelocity() {
        if (bodyOfInfluence.hasOrbit) return bodyOfInfluence.GetOrbit().GetVelocity() + velocity;
        else return velocity;
    }
    public Vector3 GetLocalVelocity() { return velocity; }
    public CelestialBody GetBodyOfInfluence() { return bodyOfInfluence; }
    //public void SetBodyOfInfluence(Transform newBoI) { bodyOfInfluence = newBoI; }
    public Orbit NewOrbit(CelestialBody newBoI) {
        Orbit newOrbit =  new Orbit(
            GetPosition() - newBoI.transform.position, // position local to newBoI
            (newBoI.hasOrbit) ? GetVelocity() - newBoI.GetOrbit().GetVelocity() : GetVelocity(), // velocity local to newBoI
            newBoI,
            mainBodyMass
        );
        return newOrbit;
    }

    public float CalculateSphereOfInfluence() { return 0.9431f * a * Mathf.Pow(mainBodyMass / bodyOfInfluence.mass, (2 / 5)); }
    //public float GetSphereOfInfluenceOf(CelestialBody bodyOfInfluence) {
    //    return 0.9431f * a * Mathf.Pow(mainBodyMass / bodyOfInfluence.mass, (2 / 5));
    //}

    // Calculates the Orbital Elements of the Orbit given a position and velocity
    public void CalculateOrbitalElementsFromPositionVelocity(Vector3 r_, Vector3 v_) {

        float r = r_.magnitude;
        float v = v_.magnitude;

        // e = Eccentricity
        this.e_ = ((v * v - sgp / r) * r_ - Vector3.Dot(r_, v_) * v_) / sgp;
        this.e = e_.magnitude;
        orbitType = FindOrbitType(e);

        // h = angular momentum
        Vector3 h_ = Vector3.Cross(r_, v_);
        h = h_.magnitude;

        // i = Inclination
        this.i = Mathf.Acos(h_.z / h);
        // Implementing a Hack for 2D 
        if (h_.z / h >= 0) this.i = 0;
        else this.i = Mathf.PI;

        // (temp) nodevector - Node Vector - points towards the ascending node
        Vector3 nodevector = Vector3.right; // Hack for 2D

        // omega = Longtitude of Ascending Node
        this.omega = Mathf.Acos(nodevector.x / nodevector.magnitude);
        if (nodevector.y < 0) this.omega = 2 * Mathf.PI - omega;
        if (i == 0 || i == Mathf.PI) {
            omega = 0;
            nodevector = Vector3.right;
            // Hack as there is no nodevector if h.z == 0
        }

        // w = argument of periapsis
        this.w = Mathf.Acos(Vector3.Dot(nodevector, e_) / (nodevector.magnitude * e));
        if (e_.z < 0) this.w = 2 * Mathf.PI - w;
        if (i == 0 || i == Mathf.PI) {
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
            position = r_; velocity = v_;
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

            float r = 2*a * (1 / (1 + Mathf.Cos(f)));

            // Flight path angle
            float phi = f/2;
            float v = Mathf.Sqrt((2 * sgp) / r);
            float vangle = f + (Mathf.PI / 2) - phi;

            P = GameSystem.PolarToCartesian(r, f);
            V = GameSystem.PolarToCartesian(v, vangle);

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
        if (inStaticOrbit) CalculateOrbitalElementsFromPositionVelocity(position, velocity);
    }

    // Draws the line of the orbit with the given number of vertices
    public void DrawOrbitalLine(LineRenderer lineRenderer, int numberOfPoints, bool pixelSnap) {
        Vector3 position;
        Vector3 lastvalidpos = Vector3.zero;
        Vector3 lastinvalidpos = Vector3.zero;

        Vector3[] positions;
        int numPositions = 0;

        bool checkIfPositionInsideSOI(Vector3 pos) {
            //float distance = (pos - bodyOfInfluence.transform.position).magnitude;
            float distance = pos.magnitude;
            return distance <= bodyOfInfluence.GetSphereOfInfluence();
        }

        Vector3 findInterceptofSOIAndLine(Vector2 pos1, Vector2 pos2) {
            //Debug.DrawLine((Vector2)bodyOfInfluence.transform.position+pos1, (Vector2)bodyOfInfluence.transform.position + pos2);
            float r = bodyOfInfluence.GetSphereOfInfluence();
            float dx = pos2.x - pos1.x;
            float dy = pos2.y - pos1.y;
            float dr = Mathf.Sqrt(dx * dx + dy * dy);
            float D = pos1.x * pos2.y - pos2.x * pos1.y;

            float discriminant = r * r * dr * dr - D * D;
            if (discriminant < 0) return Vector3.zero;

            float sgn(float n) { return (n < 0) ? -1 : 1; }

            float x1 = (D * dy + sgn(dy) * dx * Mathf.Sqrt(discriminant)) / (dr * dr);
            float x2 = (D * dy - sgn(dy) * dx * Mathf.Sqrt(discriminant)) / (dr * dr);

            float y1 = (-D * dx + Mathf.Abs(dy) * Mathf.Sqrt(discriminant)) / (dr * dr);
            float y2 = (-D * dx - Mathf.Abs(dy) * Mathf.Sqrt(discriminant)) / (dr * dr);

            //Debug.Log(string.Format("{0} {1} {2} {3}", x1, x2, y1, y2));
            //Debug.Log(pos1 + " " + pos2);
            if ((x1 >= 0) == (pos1.x+(pos2.x - pos1.x) >= 0) && (y1 >= 0) == (pos1.y+(pos2.y - pos1.y) >= 0)) return new Vector3(x1, y1, 0f);
            else return new Vector3(x2, y2, 0f);
        }

        void posCheck() {
            if (checkIfPositionInsideSOI(position)) {
                if (numPositions == 0 && lastinvalidpos != Vector3.zero) { // Entering SOI
                    positions[numPositions++] = bodyOfInfluence.transform.position + findInterceptofSOIAndLine(lastinvalidpos, position);
                }
                positions[numPositions++] = bodyOfInfluence.transform.position + position;
                lastvalidpos = position;
            }
            else if (numPositions >= 1 && positions[numPositions - 1] == bodyOfInfluence.transform.position + lastvalidpos) { // Exiting SOI
                positions[numPositions++] = bodyOfInfluence.transform.position + findInterceptofSOIAndLine(lastvalidpos, position);
            }
            else {
                lastinvalidpos = position;
            }
        }

        if (orbitType==OrbitType.Circular||orbitType==OrbitType.Elliptical) {
            numberOfPoints *= (a/3) > 9 ? (int)(Mathf.Sqrt(a) / 3) : 1;
            lineRenderer.positionCount = numberOfPoints;
            positions = new Vector3[numberOfPoints];

            for (int i = 0; i < numberOfPoints; i++) {
                position = CalculatePositionVelocityatTime(((T / (float)numberOfPoints) * (float)i + t0 + T/2), false, true);
                //lineRenderer.SetPosition(i, position);
                posCheck();
            }

            lineRenderer.SetPositions(positions);
            lineRenderer.positionCount = numPositions;
            lineRenderer.loop = numPositions == numberOfPoints;
        } else if (orbitType == OrbitType.Hyperbolic || orbitType == OrbitType.Parabolic) {
            float timerange = 50; // seconds
            lineRenderer.loop = false;
            lineRenderer.positionCount = numberOfPoints;
            positions = new Vector3[numberOfPoints];

            float thing = 2 * timerange * timerange / (numberOfPoints * numberOfPoints);
            for (int i = 1; i < numberOfPoints/2; i++) {
                float timeslice = (timerange / (float)numberOfPoints) * (timerange / ((float)i)) - thing;
                position = CalculatePositionVelocityatTime((t0+timeslice), false, true);
                posCheck();
            }
            position = bodyOfInfluence.transform.position + CalculatePositionVelocityatTime((t0), false, true);
            lineRenderer.SetPosition(numberOfPoints/2-1, position);
            for (int i = numberOfPoints / 2; i >= 1; i--) {
                float timeslice = (timerange / (float)numberOfPoints) * (timerange / ((float)i)) - thing;
                position = CalculatePositionVelocityatTime((t0 - timeslice), false, true);
                posCheck();
            }

            if (pixelSnap) for (int i = 0; i < positions.Length; i++) positions[i] = GameSystem.VPixelSnap(positions[i]);

            lineRenderer.SetPositions(positions);
            lineRenderer.positionCount = numPositions;
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
            //float rp = -a * (e - 1);
            //float n = 2 * Mathf.Sqrt(sgp / rp*rp*rp);
            float n = Mathf.Sqrt(sgp);
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
            float D = Mathf.Sqrt(l) * Mathf.Tan(f / 2);

            float M = q * D + (D * D * D) / 6f;
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
                for (int i = 0; i < 1000; i++) {
                    var dE = (E - e * Mathf.Sin(E) - M) / (1 - e * Mathf.Cos(E));
                    E -= dE;
                    if (Mathf.Abs(dE) < 1e-6) break;
                }
            }
            return E;
        } else if (orbitType == OrbitType.Parabolic) {
            // not neccessary

            return float.NaN;

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
            //if (fcompute < -1 || fcompute > 1) {
            //    Debug.LogError(string.Format("True anomaly fcompute is fucked btw (< -1); fcompute = {0}", fcompute));
            //}
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
            // Use Barker's equation to calcluate true anomaly
            float rp = a;
            float A = (3f / 2) * Mathf.Sqrt((sgp) / (2 * rp * rp * rp)) * (currentTime - t0);
            float B = Mathf.Pow((A + Mathf.Sqrt(A * A + 1)), (1 / 3f));
            float newf = 2 * Mathf.Atan(B - (1 / B));
            return newf;
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
        }
        else if (orbitType == OrbitType.Parabolic) {
            // q = a, cause rp = a = p/2 = q;
            float tempf = 2 * Mathf.Atan(E/Mathf.Sqrt(2*a));
            return tempf;
        }
        else if (orbitType == OrbitType.Hyperbolic) {
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

    public OrbitType FindOrbitType(float e) {
        if (e == 0) return OrbitType.Circular;
        else if (e > 0 && e < 1) return OrbitType.Elliptical;
        else if (e == 1) return OrbitType.Parabolic;
        else if (e > 1) return OrbitType.Hyperbolic;
        else return 0;
    }

    public override string ToString() {
        string outtext = "";
        outtext += string.Format("BoI: {0}\n", bodyOfInfluence.name);
        outtext += string.Format("Position, Velocity: {0}, {1}\n", GetLocalPosition(), GetLocalVelocity());
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
        outtext += string.Format("Periapsis (m): {0}\n", pe);
        outtext += string.Format("Apoapsis (m): {0}\n", ap);

        return outtext;
    }
}
