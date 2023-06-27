using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RocketScript : MonoBehaviour {

    public UIController UIController;
    private LineRenderer lineRenderer;
    public Rigidbody2D rigidbody2D;

    public CelestialBody bodyOfInfluence;

    public OrbitDouble rocketOrbit;

    public float mass = 1f;

    public float thrustPower = 1f;

    [Range(0f, 10f)]
    public float rotateSpeed;

    [Range(0, 200)]
    public int orbitalLineNumberOfPoints;

    private bool drawOrbitalLine = true;

    [Header("Sprites")]
    [Space(10)]
    public Sprite normalSprite;
    public Sprite boostSprite;

    public Vector3 startPositionVector = new Vector3(6f, 0f);
    public Vector3 startVelocityVector = new Vector3(-1f, 15f);

    private bool isThrust = false;

    // Start is called before the first frame update
    void Start() {
        rocketOrbit = new OrbitDouble(startPositionVector, startVelocityVector, bodyOfInfluence, mass);
        lineRenderer = GetComponent<LineRenderer>();
        rigidbody2D = GetComponent<Rigidbody2D>();
    }

    // Update is called once per frame
    public void OrbitalUpdate() {

        GameSystem.timeElapsed += GameSystem.DeltaTime();

        if (Input.GetKey(KeyCode.Q)) GameSystem.Rotate(transform, rotateSpeed * Time.deltaTime);
        if (Input.GetKey(KeyCode.E)) GameSystem.Rotate(transform, -rotateSpeed * Time.deltaTime);

        if (Input.GetKeyDown(KeyCode.B)) {
            if (rocketOrbit.inStaticOrbit) rocketOrbit.ChangeOrbitToDynamic(ref rigidbody2D);
            else rocketOrbit.ChangeOrbitToStatic();
        }

        if (rocketOrbit.inStaticOrbit) {
            
            rocketOrbit.CalculatePositionVelocityatTime(GameSystem.CurrentTime());
            transform.position = rocketOrbit.GetPosition();

            rocketOrbit.CalculateExtraVariables();

        }
        // Draw orbital line
        if (drawOrbitalLine) {
            lineRenderer.enabled = true;
            GameSystem.SetLineWidth(lineRenderer, GameSystem.pixelUnit * GameSystem.screenScale);
            rocketOrbit.DrawOrbitalLine(lineRenderer, orbitalLineNumberOfPoints, false);
        }
        else lineRenderer.enabled = false;

        if (GameSystem.DEBUG) UIController.SetText(rocketOrbit.ToString() + "\nTimeScale: " + GameSystem.timeScale + "\ndrawOrbitalLine: " + drawOrbitalLine);

        if (isThrust) SetSprite(boostSprite);
        else SetSprite(normalSprite);
        //SetSpriteSize();

        // Check which body of influence it is within the SOI of
        CelestialBody currentBOI = CheckSOIInterceptionOfChildren(bodyOfInfluence.GetRootBody());
        if (!currentBOI.Equals(bodyOfInfluence)) {
            bodyOfInfluence = currentBOI;
            rocketOrbit = rocketOrbit.NewOrbit(currentBOI);
            Debug.Log("new BOI = " + currentBOI);
            GameSystem.timeScale = 1f;
        }

    }
    public void FixedUpdate()
    {
        if (!rocketOrbit.inStaticOrbit) {
            rocketOrbit.SetPosition(rigidbody2D.position);
            rocketOrbit.CalculateOrbitalElementsFromPositionVelocity(rocketOrbit.GetLocalPosition(), rocketOrbit.GetLocalVelocity());
        }

        isThrust = false;
        if (Input.GetKey(KeyCode.W)) {
            rocketOrbit.AddForce(thrustPower * GameSystem.AngleToDirection(transform.rotation.eulerAngles.z), GameSystem.FixedDeltaTime(), rocketOrbit.inStaticOrbit);
            isThrust = true;
        }

        if (!rocketOrbit.inStaticOrbit) {
            rocketOrbit.AddForce(Vector3.zero, GameSystem.FixedDeltaTime(), true); // Only adds gravity
            rigidbody2D.MovePosition(rocketOrbit.GetPosition());
        }
    }
    public void OnCollisionEnter2D(Collision2D collision)
    {
        //rocketOrbit.SetVelocity(Vector2.zero);
        drawOrbitalLine = false;
        Debug.Log("drawOrbitalLine = false");
    }
    public void OnCollisionExit2D(Collision2D collision) {
        drawOrbitalLine = true;
        Debug.Log("drawOrbitalLine = true");
    }

    void SetSprite(Sprite newsprite) {
        GetComponent<SpriteRenderer>().sprite = newsprite;
    }
    void SetSpriteSize() {
        GetComponent<SpriteRenderer>().size = Vector2.one * GameSystem.screenScale;
    }

    private CelestialBody CheckSOIInterceptionOfChildren(CelestialBody rootBody) {
        CelestialBody highest = rootBody;
        int highestnum = 0;
        foreach (Transform bodytransform in rootBody.transform) {
            CelestialBody body = bodytransform.GetComponent<CelestialBody>();
            int heirarchy = body.GetHeirarchyLevel();
            if (body.GetHeirarchyLevel() > highestnum) {
                float distance = (transform.position - bodytransform.position).magnitude;
                if (distance <= body.GetSphereOfInfluence()) {
                    if (!body.Equals(bodyOfInfluence) && GameSystem.DEBUG) Debug.Log("within radius of " + body);
                    highest = body;
                    highestnum = heirarchy;
                    if (bodytransform.childCount > 0) {
                        CelestialBody highestchild = CheckSOIInterceptionOfChildren(body);
                        if (highestchild.GetHeirarchyLevel() > heirarchy) return highestchild;
                    }
                }
            }
        }
        return highest;
    }
}
