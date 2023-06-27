using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UIElements;

public class OrbitDouble {

    public enum OrbitType { Circular, Elliptical, Parabolic, Hyperbolic }

    public static double gravconst = 6.6725985E-11f; // fundamental universal constant

    private CelestialBody bodyOfInfluence;

    private Vector3 position;
    private Vector3 velocity;
    private double mainBodyMass;
    private double bodyOfInfluenceMass;

    public bool inStaticOrbit = true;

    public OrbitType orbitType;

    // Standard gravitational parameter - Mu
    private double sgp;

    // Specific Angular Momentum - h (conserved in closed system)
    private double h;

    // == Orbital Elements of the Object's orbit == 

    // Eccentricity Vector
    private Vector3 e_; 

    // Eccentricity
    private double e;

    // Semi-Major Axis (m)
    private double a;

    // Argument of periapsis (rad)
    private double w;

    // Orbital Period
    private double T;

    // Epoch - Always set to time at periapsis
    private double t0;

    // Inclination
    private double i;

    // Longtitude of Right Ascending Node (RAAN)
    private double omega;

    // True Anomaly (rad)
    private double f;

    // == Extra variables ==

    private Vector3 pe_;
    private Vector3 ap_;
    private double pe;
    private double ap;


    // Initialises the state variables in the orbit, given BoI, and calculates orbital elements
    public OrbitDouble(Vector3 position, Vector3 velocity, CelestialBody bodyOfInfluence, double mainBodyMass) {
        this.mainBodyMass = mainBodyMass;
        this.bodyOfInfluenceMass = bodyOfInfluence.mass;
        this.bodyOfInfluence = bodyOfInfluence;

        this.sgp = bodyOfInfluenceMass * gravconst;

        CalculateOrbitalElementsFromPositionVelocity(position, velocity);
    }

    // Constructor that gives us the 6 Orbital Elements + BOI Information
    public OrbitDouble(double eccentricity, double semiMajorAxis, double argumentOfPeriapsis,
        double inclination, double longtitudeOfAscendingNode, double trueAnomaly,
        CelestialBody bodyOfInfluence, double mainBodyMass
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

        double n = CalculateMeanMotion(orbitType);

        // T = Orbital Period
        if (orbitType == OrbitType.Circular || orbitType == OrbitType.Elliptical) {
            this.T = (2 * System.Math.PI) / n;
        }
        else this.T = float.PositiveInfinity;

        if (f != 0) this.t0 = GameSystem.CurrentTime() - CalculateMeanAnomaly(orbitType) / n;
        else this.t0 = 0;
        //Debug.Log(M);

        CalculatePositionVelocityatTime(GameSystem.CurrentTime());

        //CalculateOrbitalElementsFromPositionVelocity(position, velocity);
    }

    public Vector3 GetPosition() { return position + ((bodyOfInfluence.hasOrbit) ? bodyOfInfluence.GetOrbit().GetPosition() : bodyOfInfluence.transform.position); }
    public Vector3 GetLocalPosition() { return position; }
    public Vector3 GetVelocity() {
        if (bodyOfInfluence.hasOrbit) return bodyOfInfluence.GetOrbit().GetVelocity() + velocity;
        else return velocity;
    }
    public Vector3 GetLocalVelocity() { return velocity; }
    public CelestialBody GetBodyOfInfluence() { return bodyOfInfluence; }

    public void SetPosition(Vector3 position) { this.position = position; }
    public void SetVelocity(Vector3 velocity) { this.velocity = velocity; }
    //public void SetBodyOfInfluence(Transform newBoI) { bodyOfInfluence = newBoI; }

    public OrbitDouble NewOrbit(CelestialBody newBoI) {
        OrbitDouble newOrbit =  new OrbitDouble(
            (newBoI.hasOrbit) ? GetPosition() - newBoI.GetOrbit().GetPosition() : GetPosition(), // position local to newBoI
            (newBoI.hasOrbit) ? GetVelocity() - newBoI.GetOrbit().GetVelocity() : GetVelocity(), // velocity local to newBoI
            newBoI,
            mainBodyMass
        );
        return newOrbit;
    }

    // Calculates the Orbital Elements of the Orbit given a position and velocity
    public void CalculateOrbitalElementsFromPositionVelocity(Vector3 r_, Vector3 v_) {

        double r = r_.magnitude;
        double v = v_.magnitude;

        // e = Eccentricity
        this.e_ = ((float)(v * v - sgp / r) * r_ - Vector3.Dot(r_, v_) * v_) / (float)sgp;
        this.e = e_.magnitude;
        orbitType = FindOrbitType(e);

        // h = angular momentum
        Vector3 h_ = Vector3.Cross(r_, v_);
        h = h_.magnitude;

        // i = Inclination
        this.i = System.Math.Acos(h_.z / h);
        // Implementing a Hack for 2D 
        if (h_.z / h >= 0) this.i = 0;
        else this.i = System.Math.PI;

        // (temp) nodevector - Node Vector - points towards the ascending node
        Vector3 nodevector = Vector3.right; // Hack for 2D

        // omega = Longtitude of Ascending Node
        this.omega = System.Math.Acos(nodevector.x / nodevector.magnitude);
        if (nodevector.y < 0) this.omega = 2 * System.Math.PI - omega;
        if (i == 0 || i == System.Math.PI) {
            omega = 0;
            nodevector = Vector3.right;
            // Hack as there is no nodevector if h.z == 0
        }

        // w = argument of periapsis
        this.w = System.Math.Acos(Vector3.Dot(nodevector, e_) / (nodevector.magnitude * e));
        if (e_.z < 0) this.w = 2 * System.Math.PI - w;
        if (i == 0 || i == System.Math.PI) {
            this.w = System.Math.Atan2(e_.y, e_.x);
            if (Vector3.Cross(r_, v_).z < 0) this.w = 2 * System.Math.PI - w;
            // Hack as there is no nodevector if h.z == 0
        }

        // a = semi-major axis
        this.a = 1 / (2 / r - (v * v) / sgp);

        this.f = CalculateTrueAnomaly(r_, v_, orbitType);

        // Find M at this point. Then use that to find t0
        // M at t0 = 0. Thus the difference in t is M/n

        double M = CalculateMeanAnomaly(orbitType);

        double n = CalculateMeanMotion(orbitType);

        // T = Orbital Period
        if (orbitType == OrbitType.Circular || orbitType == OrbitType.Elliptical) {
            this.T = (2 * System.Math.PI) / n;
        }
        else this.T = float.PositiveInfinity;

        t0 = GameSystem.CurrentTime() - (M / n); // Sets the epoch to the time of periapsis


        if (GameSystem.LOG_ELEMENTS && inStaticOrbit) {
            position = r_; velocity = v_;
            Debug.Log(ToString());
        }


    }

    // Calculates the position and velocity of mainBody given a time
    public Vector3 CalculatePositionVelocityatTime(double inputTime, bool setvariables=true, bool evenmotion=false) {

        double n = CalculateMeanMotion(orbitType);
        double f = CalculateTrueAnomalyFromTime(inputTime, orbitType);
        //if (orbitType == OrbitType.Hyperbolic && setvariables) Debug.Log("f at " + Time.frameCount + " is " + f);

        Vector2 P = Vector2.zero, V = Vector2.zero;


        if (orbitType == OrbitType.Circular || orbitType == OrbitType.Elliptical) {
            double E = CalculateEccentricAnomaly(inputTime, orbitType, evenmotion);
            // Calculate X, Y and Z values
            P.x = (float)(a * (System.Math.Cos(E) - e));
            P.y = (float)(a * System.Math.Sin(E) * System.Math.Sqrt(1 - e * e));

            // Calculate vX, vY and vZ values
            V.x = (float)(-a * System.Math.Sin(E) * n / (1 - e * System.Math.Cos(E)));
            V.y = (float)(a * System.Math.Cos(E) * System.Math.Sqrt(1 - e * e) * n / (1 - e * System.Math.Cos(E)));
        } else if (orbitType == OrbitType.Parabolic) {
            double r = 2*a * (1 / (1 + System.Math.Cos(f)));

            // Flight path angle
            double phi = f/2;
            double v = System.Math.Sqrt((2 * sgp) / r);
            double vangle = f + (System.Math.PI / 2) - phi;

            P = GameSystem.PolarToCartesian((float)r, (float)f);
            V = GameSystem.PolarToCartesian((float)v, (float)vangle);

        } else if (orbitType == OrbitType.Hyperbolic) {
            double l = a * (e * e - 1);
            double r = -(l / (1 + e * System.Math.Cos(f))); // idk why tf but it has to be negative

            // Flight path angle
            double phi = System.Math.Atan((e * System.Math.Sin(f)) / (1 + e * System.Math.Cos(f)));
            double v = System.Math.Sqrt(sgp * ((2 / r) - (1 / a)));
            double vangle = f + (System.Math.PI / 2) - phi;

            P = GameSystem.PolarToCartesian((float)r, (float)f);
            V = GameSystem.PolarToCartesian((float)v, (float)vangle);
        }

        Vector3 r_ = RotateFromOrbitalPlaneTo3D(P.x, P.y);
        Vector3 v_ = RotateFromOrbitalPlaneTo3D(V.x, V.y);

        if (setvariables) {
            this.f = f;
            position = r_;
            velocity = v_;
            return Vector3.zero;
        }
        else return r_;

    }

    // Simulate the force of gravity for the current physics interval
    public Vector3 SimulateGravity(float timeStep) {
        Vector3 difference = Vector3.zero - position;
        float r = difference.magnitude;
        Vector3 gravacc = (float)(timeStep * (sgp / (r * r))) * difference.normalized;
        return gravacc;
    }
    
    public void ChangeOrbitToDynamic(ref Rigidbody2D rigidbody2D) {
        rigidbody2D.position = position;
        rigidbody2D.velocity = velocity;
        inStaticOrbit = false;
    }
    public void ChangeOrbitToStatic() {
        inStaticOrbit = true;
    }

    // Adds the acceleration of a given force and recalculates orbital elements - also adds step in position
    public void AddForce(Vector3 force, float timeStep, bool doGravity) {
        Vector3 thrustacc = timeStep * force / (float)mainBodyMass;
        // a = (GM)/r^2

        if (doGravity) velocity += SimulateGravity(timeStep);
        velocity += thrustacc;
        position += velocity * timeStep; // Adding acceleration of forces
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
            //double distance = (pos - bodyOfInfluence.transform.position).magnitude;
            double distance = pos.magnitude;
            return distance <= bodyOfInfluence.GetSphereOfInfluence();
        }

        Vector3 findInterceptofSOIAndLine(Vector2 pos1, Vector2 pos2) {
            //Debug.DrawLine((Vector2)bodyOfInfluence.transform.position+pos1, (Vector2)bodyOfInfluence.transform.position + pos2);
            double r = bodyOfInfluence.GetSphereOfInfluence();
            double dx = pos2.x - pos1.x;
            double dy = pos2.y - pos1.y;
            double dr = System.Math.Sqrt(dx * dx + dy * dy);
            double D = pos1.x * pos2.y - pos2.x * pos1.y;

            double discriminant = r * r * dr * dr - D * D;
            if (discriminant < 0) return Vector3.zero;

            double sgn(double n) { return (n < 0) ? -1 : 1; }

            double x1 = (D * dy + sgn(dy) * dx * System.Math.Sqrt(discriminant)) / (dr * dr);
            double x2 = (D * dy - sgn(dy) * dx * System.Math.Sqrt(discriminant)) / (dr * dr);

            double y1 = (-D * dx + System.Math.Abs(dy) * System.Math.Sqrt(discriminant)) / (dr * dr);
            double y2 = (-D * dx - System.Math.Abs(dy) * System.Math.Sqrt(discriminant)) / (dr * dr);

            //Debug.Log(string.Format("{0} {1} {2} {3}", x1, x2, y1, y2));
            //Debug.Log(pos1 + " " + pos2);
            if ((x1 >= 0) == (pos1.x+(pos2.x - pos1.x) >= 0) && (y1 >= 0) == (pos1.y+(pos2.y - pos1.y) >= 0)) return new Vector3((float)x1, (float)y1, 0f);
            else return new Vector3((float)x2, (float)y2, 0f);
        }

        void posCheck() {
            if (checkIfPositionInsideSOI(position)) {
                if (numPositions == 0 && lastinvalidpos != Vector3.zero) { // Entering SOI
                    positions[numPositions++] = bodyOfInfluence.transform.position + findInterceptofSOIAndLine(lastinvalidpos, position);
                    //if (!checkIfPositionInsideSOI(positions[numPositions])) Debug.LogError("Position on circle outside SOI");
                }
                positions[numPositions++] = bodyOfInfluence.transform.position + position;
                lastvalidpos = position;
            }
            else if (numPositions >= 1 && positions[numPositions - 1] == bodyOfInfluence.transform.position + lastvalidpos) { // Exiting SOI
                positions[numPositions++] = bodyOfInfluence.transform.position + findInterceptofSOIAndLine(lastvalidpos, position);
                //if (!checkIfPositionInsideSOI(positions[numPositions])) Debug.LogError("Position on circle outside SOI");
            }
            else {
                lastinvalidpos = position;
            }
        }

        if (orbitType==OrbitType.Circular||orbitType==OrbitType.Elliptical) {
            numberOfPoints *= (a/3) > 9 ? (int)(System.Math.Sqrt(a) / 3) : 1;
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
            double timerange = 5000000; // seconds
            lineRenderer.loop = false;
            lineRenderer.positionCount = numberOfPoints;
            positions = new Vector3[numberOfPoints];

            double thing = 2 * timerange * timerange / (numberOfPoints * numberOfPoints);
            for (int i = 1; i < numberOfPoints/2; i++) {
                double timeslice = (timerange / (float)numberOfPoints) * (timerange / ((float)i)) - thing;
                position = CalculatePositionVelocityatTime((t0+timeslice), false, true);
                posCheck();
            }
            position = bodyOfInfluence.transform.position + CalculatePositionVelocityatTime((t0), false, true);
            lineRenderer.SetPosition(numberOfPoints/2-1, position);
            for (int i = numberOfPoints / 2; i >= 1; i--) {
                double timeslice = (timerange / (float)numberOfPoints) * (timerange / ((float)i)) - thing;
                position = CalculatePositionVelocityatTime((t0 - timeslice), false, true);
                posCheck();
            }

            if (pixelSnap) for (int i = 0; i < positions.Length; i++) positions[i] = GameSystem.VPixelSnap(positions[i]);

            lineRenderer.SetPositions(positions);
            lineRenderer.positionCount = numPositions;
        }
    }

    Vector3 RotateFromOrbitalPlaneTo3D(double P, double Q) {

        // rotate by argument of periapsis
        double x = System.Math.Cos(w) * P - System.Math.Sin(w) * Q;
        double y = System.Math.Sin(w) * P + System.Math.Cos(w) * Q;
        // rotate by inclination
        double z = System.Math.Sin(i) * y;
        y = System.Math.Cos(i) * y;
        // rotate by longitude of ascending node
        double xtemp = x;
        x = System.Math.Cos(omega) * xtemp - System.Math.Sin(omega) * y;
        y = System.Math.Sin(omega) * xtemp + System.Math.Cos(omega) * y;

        return new Vector3((float)x, (float)y, (float)z);
    }

    double CalculateMeanMotion(OrbitType orbitType) {
        if (orbitType==OrbitType.Circular||orbitType==OrbitType.Elliptical) {
            double n = System.Math.Sqrt(sgp / (a * a * a));
            return n;
        } else if (orbitType==OrbitType.Parabolic) {
            //double rp = -a * (e - 1);
            //double n = 2 * System.Math.Sqrt(sgp / rp*rp*rp);
            double n = System.Math.Sqrt(sgp);
            return n;
        }
        else if (orbitType==OrbitType.Hyperbolic) {
            double n = System.Math.Sqrt(sgp / -(a * a * a));
            return n;
        }
        return 0;
    }

    double CalculateMeanAnomaly(OrbitType orbitType) {
        // This assumes f is given

        if (orbitType==OrbitType.Circular||orbitType==OrbitType.Elliptical) {
            double M = System.Math.Atan2(-System.Math.Sqrt(1 - e * e) * System.Math.Sin(f), -e - System.Math.Cos(f)) + System.Math.PI - e * ((System.Math.Sqrt(1 - e * e) * System.Math.Sin(f)) / (1 + e * System.Math.Cos(f)));
            return M;
        } else if (orbitType == OrbitType.Parabolic) {
            double l = a * (e * e - 1);
            double q = l / 2;
            double D = System.Math.Sqrt(l) * System.Math.Tan(f / 2);

            double M = q * D + (D * D * D) / 6f;
            return M;
        } else if (orbitType==OrbitType.Hyperbolic) {
            double first = System.Math.Tan(f / 2) / System.Math.Sqrt(((e + 1) / (e - 1)));
            double F = 2 * System.Math.Atanh(first);

            double M = (e * System.Math.Sinh(F)) - F;
            return M;
        }
        return 0;
    }

    double CalculateEccentricAnomaly(double currentTime, OrbitType orbitType, bool evenmotion=false) {

        double timeSincePeriapsis = currentTime - t0;

        if (orbitType==OrbitType.Circular) {
            double n = CalculateMeanMotion(orbitType);
            double M = n * (timeSincePeriapsis); // Mean Anomaly (rad)
            return M;
        }
        else if (orbitType == OrbitType.Elliptical) {
            // Using Newton's approximation method
            double n = CalculateMeanMotion(orbitType);
            double M = n * (timeSincePeriapsis); // Mean Anomaly (rad)

            double E = M; // Eccentric Anomaly (rad) - E
            if (!evenmotion) {
                for (int i = 0; i < 1000; i++) {
                    var dE = (E - e * System.Math.Sin(E) - M) / (1 - e * System.Math.Cos(E));
                    E -= dE;
                    if (System.Math.Abs(dE) < 1e-6) break;
                }
            }
            return E;
        } else if (orbitType == OrbitType.Parabolic) {
            // not neccessary

            return float.NaN;

        } else if (orbitType == OrbitType.Hyperbolic) {
            // Using Newton's approximation method
            double n = CalculateMeanMotion(orbitType);
            double M = n * (timeSincePeriapsis); // Mean Anomaly (rad)

            double F = M; // Eccentric Anomaly (rad) - E
            for (int i = 0; i < 1000; i++) {
                var dF = (e * System.Math.Sinh(F) - F - M) / (e * System.Math.Cosh(F)-1);
                F -= dF;
                if (System.Math.Abs(dF) < 1e-6) break;
            }
            return F;
        }

        return 0;

    }

    double CalculateTrueAnomaly(Vector3 r_, Vector3 v_, OrbitType orbitType) {
        // We use argument of latitude if eccentricity is 0 (circular orbit)

        if (orbitType == OrbitType.Circular) { // Circular Orbit
            Vector3 nodevector = Vector3.right; // Hack for 2D

            double u = System.Math.Acos((Vector3.Dot(nodevector, r_) / (nodevector.magnitude * r_.magnitude)));
            if (r_.z<0) u = 2 * System.Math.PI - u;

            return u;
        }else {
            double fcompute = Vector3.Dot(e_, r_) / (e * r_.magnitude);
            //if (fcompute < -1 || fcompute > 1) {
            //    Debug.LogError(string.Format("True anomaly fcompute is fucked btw (< -1); fcompute = {0}", fcompute));
            //}
            double tempf = System.Math.Acos(fcompute);
            if (Vector3.Dot(r_, v_) < 0) tempf = 2 * System.Math.PI - tempf;

            //double h = System.Math.Sqrt(sgp * (a * (e * e - 1)));
            Vector3 h_ = Vector3.Cross(r_, v_);
            double h = h_.magnitude;

            double q = Vector3.Dot(r_, v_);
            double fy = (h * q) / (r_.magnitude * sgp);
            double fx = (h * h) / (r_.magnitude * sgp) - 1;
            double tempf2 = System.Math.Atan2(fy, fx);

            return tempf2;
        }

    }

    double CalculateTrueAnomalyFromTime(double currentTime, OrbitType orbitType) {
        // We use argument of latitude if eccentricity is 0 (circular orbit)

        if (orbitType == OrbitType.Circular) { // Circular Orbit
            double E = CalculateEccentricAnomaly(currentTime, orbitType, false);
            return E;
        }
        else if (orbitType == OrbitType.Elliptical) { // Eliptical Orbit
            double E = CalculateEccentricAnomaly(currentTime, orbitType, false);
            return CalculateTrueAnomalyFromE(E, orbitType);
        }
        else if (orbitType == OrbitType.Parabolic) { // Parabolic Orbit
            // Use Barker's equation to calcluate true anomaly
            double rp = a;
            double A = (3f / 2) * System.Math.Sqrt((sgp) / (2 * rp * rp * rp)) * (currentTime - t0);
            double B = System.Math.Pow((A + System.Math.Sqrt(A * A + 1)), (1 / 3f));
            double newf = 2 * System.Math.Atan(B - (1 / B));
            return newf;
        }
        else if (orbitType == OrbitType.Hyperbolic) { // Hyperbolic Orbit
            double F = CalculateEccentricAnomaly(currentTime, orbitType, false);
            //if (currentTime==GameSystem.CurrentTime()) Debug.Log("F at " + Time.frameCount + " is " + F);
            return CalculateTrueAnomalyFromE(F, orbitType);
        }

        return 0;

    }

    double CalculateTrueAnomalyFromE(double E, OrbitType orbitType) {
        if (orbitType==OrbitType.Elliptical) {
            double B = e / (1 + System.Math.Sqrt(1 - e * e));
            double tempf = E + 2 * System.Math.Atan((B * System.Math.Sin(E)) / (1 - B * System.Math.Cos(E)));
            return tempf;
        }
        else if (orbitType == OrbitType.Parabolic) {
            // q = a, cause rp = a = p/2 = q;
            double tempf = 2 * System.Math.Atan(E/System.Math.Sqrt(2*a));
            return tempf;
        }
        else if (orbitType == OrbitType.Hyperbolic) {
            double tempf = 2*System.Math.Atan(System.Math.Sqrt((e + 1) / (e - 1)) * System.Math.Tanh(E / 2));
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

    public OrbitType FindOrbitType(double e) {
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
        outtext += string.Format("inStaticOrbit: {0}\n", inStaticOrbit);

        return outtext;
    }
}
