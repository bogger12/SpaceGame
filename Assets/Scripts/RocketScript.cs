using System.Collections;
using System.Collections.Generic;
using UnityEngine;
//using static UnityEditor.PlayerSettings;

public class RocketScript : MonoBehaviour {

    public UIController UIController;
    public LineRenderer lineRenderer;

    public CelestialBody bodyOfInfluence;

    public Orbit rocketOrbit;

    private float mass = 1f;

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
        //rigidbody2D = GetComponent<Rigidbody2D>();
    }
    
    // Update is called once per frame
    private void Update() {
        isThrust = false;
        if (Input.GetKey(KeyCode.W)) {
            rocketOrbit.AddForce(GameSystem.AngleToDirection(transform.rotation.eulerAngles.z));
            isThrust = true;
        }

        //if (Input.GetKey(KeyCode.Q)) rigidbody2D.MoveRotation(rigidbody2D.rotation + Mathf.Rad2Deg * rotateSpeed * Time.deltaTime);
        //if (Input.GetKey(KeyCode.E)) rigidbody2D.MoveRotation(rigidbody2D.rotation + Mathf.Rad2Deg * -rotateSpeed * Time.deltaTime);
        if (Input.GetKey(KeyCode.Q)) GameSystem.Rotate(transform, rotateSpeed * Time.deltaTime);
        if (Input.GetKey(KeyCode.E)) GameSystem.Rotate(transform, -rotateSpeed * Time.deltaTime);

        if (!rocketOrbit.landed) {
            rocketOrbit.CalculatePositionVelocityatTime(Time.time);
            //transform.position = GameSystem.VPixelSnap(rocketOrbit.GetPosition());
            transform.position = rocketOrbit.GetPosition();
            //rigidbody2D.MovePosition(rocketOrbit.GetPosition());
            //Debug.Log("rocket at " + transform.position);
            rocketOrbit.CalculateExtraVariables();

            // Draw orbital line
            GameSystem.SetLineWidth(lineRenderer, GameSystem.pixelUnit * GameSystem.screenScale);
            rocketOrbit.DrawOrbitalLine(lineRenderer, orbitalLineNumberOfPoints, false);

            //Debug.Log(rocketOrbit);
            if (GameSystem.DEBUG) UIController.SetText(rocketOrbit.ToString());

            CelestialBody currentBOI = CheckSOIInterceptionOfChildren(bodyOfInfluence.GetRootBody());
            if (!currentBOI.Equals(bodyOfInfluence)) {
                bodyOfInfluence = currentBOI;
                rocketOrbit = rocketOrbit.NewOrbit(currentBOI);
                Debug.Log("new BOI = " + currentBOI);
            }
        }
        else {
            transform.position = rocketOrbit.GetPosition();
        }
        
        if (isThrust) SetSprite(boostSprite);
        else SetSprite(normalSprite);

    }

    void SetSprite(Sprite newsprite) {
        GetComponent<SpriteRenderer>().sprite = newsprite;
    }

    //private void OnTriggerEnter2D(Collider2D collision)
    //{
    //    //CircleCollider2D collider = collision.GetComponent<CircleCollider2D>();
    //    Debug.Log(collision);
    //    Debug.Log("CollisionEnter of " + collision.gameObject + " at " + transform.position + ", collisionobj at " + collision.gameObject.transform.position);
    //    if (!collision.Equals(bodyOfInfluence.sphereOfInfluence)) {
    //        CelestialBody colliderBodyOfInfluence = collision.gameObject.GetComponent<CelestialBody>();
    //        Debug.Log("newCollisionEnter");
    //        if (collision.Equals(colliderBodyOfInfluence.sphereOfInfluence)) {
    //            bodyOfInfluence = colliderBodyOfInfluence;
    //            Debug.Log("newBodyOfInfluence = " + colliderBodyOfInfluence);
    //            rocketOrbit = rocketOrbit.NewOrbit(colliderBodyOfInfluence);
    //        }
    //    }
    //}

    private CelestialBody CheckSOIInterceptionOfChildren(CelestialBody rootBody) {
        CelestialBody highest = rootBody;
        int highestnum = 0;
        foreach (Transform bodytransform in rootBody.transform) {
            CelestialBody body = bodytransform.GetComponent<CelestialBody>();
            int heirarchy = body.GetHeirarchyLevel();
            if (body.GetHeirarchyLevel() > highestnum) {
                float distance = (transform.position - bodytransform.position).magnitude;
                if (distance <= body.radiusSOI) {
                    if (!body.Equals(bodyOfInfluence)) Debug.Log("within radius of " + body);
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
