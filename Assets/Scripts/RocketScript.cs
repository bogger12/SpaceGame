using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RocketScript : MonoBehaviour {

    public UIController UIController;
    public LineRenderer lineRenderer;

    public CelestialBody bodyOfInfluence;

    public Orbit rocketOrbit;

    public float mass = 1f;

    [Range(0f, 10f)]
    public float rotateSpeed;

    [Range(0, 200)]
    public int orbitalLineNumberOfPoints;

    [Header("Sprites")]
    [Space(10)]
    public Sprite normalSprite;
    public Sprite boostSprite;

    public Vector3 startPositionVector = new Vector3(6f, 0f);
    public Vector3 startVelocityVector = new Vector3(-1f, 15f);

    private bool isThrust = false;

    // Start is called before the first frame update
    void Start() {
        rocketOrbit = new Orbit(startPositionVector, startVelocityVector, bodyOfInfluence, mass);
    }
    
    // Update is called once per frame
    private void Update() {
        isThrust = false;
        if (Input.GetKey(KeyCode.W)) {
            rocketOrbit.AddForce(GameSystem.AngleToDirection(transform.rotation.eulerAngles.z));
            isThrust = true;
        }

        if (Input.GetKey(KeyCode.Q)) GameSystem.Rotate(transform, rotateSpeed * Time.deltaTime);
        if (Input.GetKey(KeyCode.E)) GameSystem.Rotate(transform, -rotateSpeed * Time.deltaTime);

        if (rocketOrbit.inStaticOrbit) {
            rocketOrbit.CalculatePositionVelocityatTime(Time.time);
            transform.position = rocketOrbit.GetPosition();

            rocketOrbit.CalculateExtraVariables();

            // Draw orbital line
            GameSystem.SetLineWidth(lineRenderer, GameSystem.pixelUnit * GameSystem.screenScale);
            rocketOrbit.DrawOrbitalLine(lineRenderer, orbitalLineNumberOfPoints, false);

            if (GameSystem.DEBUG) UIController.SetText(rocketOrbit.ToString());

            // Check which body of influence it is within the SOI of
            CelestialBody currentBOI = CheckSOIInterceptionOfChildren(bodyOfInfluence.GetRootBody());
            if (!currentBOI.Equals(bodyOfInfluence)) {
                bodyOfInfluence = currentBOI;
                rocketOrbit = rocketOrbit.NewOrbit(currentBOI);
                Debug.Log("new BOI = " + currentBOI);
            }
        }
        else {
            
            //transform.position = rocketOrbit.GetPosition();
        }

        if (isThrust) SetSprite(boostSprite);
        else SetSprite(normalSprite);

    }

    void SetSprite(Sprite newsprite) {
        GetComponent<SpriteRenderer>().sprite = newsprite;
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
