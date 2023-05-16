using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RocketScript : MonoBehaviour {

    public UIController UIController;
    public LineRenderer lineRenderer;
    public GameObject bodyOfInfluence;

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

    public Vector3 startVector = new Vector3(-1f, 15f);


    private bool isThrust = false;

    // Start is called before the first frame update
    void Start() {
        rocketOrbit = new Orbit(transform.position, startVector, mass, 5.9722E12f);
    }
    
    // Update is called once per frame
    private void Update() {
        isThrust = false;
        if (Input.GetKey(KeyCode.W)) {
            rocketOrbit.AddForce(GetDirection());
            isThrust = true;
        }

        if (Input.GetKey(KeyCode.Q)) GameSystem.Rotate(transform, rotateSpeed * Time.deltaTime);
        if (Input.GetKey(KeyCode.E)) GameSystem.Rotate(transform, -rotateSpeed * Time.deltaTime);

        rocketOrbit.CalculatePositionVelocityatTime(Time.time);
        transform.position = rocketOrbit.position;

        rocketOrbit.DrawOrbitalLine(lineRenderer, orbitalLineNumberOfPoints);

        UIController.SetText(rocketOrbit.ToString());

        rocketOrbit.CalculateExtraVariables();

        if (isThrust) SetSprite(boostSprite);
        else SetSprite(normalSprite);
    }

    Vector3 GetDirection() {
        float dirangle = Mathf.Deg2Rad*(transform.rotation.eulerAngles.z+90);
        //Debug.Log(dirangle);
        float x = Mathf.Cos(dirangle);
        float y = Mathf.Sin(dirangle);

        return new Vector3(x, y);
    }

    void SetSprite(Sprite newsprite) {
        GetComponent<SpriteRenderer>().sprite = newsprite;
    }

}
