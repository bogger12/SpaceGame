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



    // Start is called before the first frame update
    void Start() {
        rocketOrbit = new Orbit(transform.position, new Vector3(-1f, 15f), mass, 5.9722E12f);
    }
    
    // Update is called once per frame
    private void Update() {
        if (Input.GetKey(KeyCode.W)) {
            rocketOrbit.AddForce(GetDirection());
        }
        if (Input.GetKey(KeyCode.S)) {
            rocketOrbit.AddForce(-GetDirection());
        }

        if (Input.GetKey(KeyCode.Q)) Rotate(rotateSpeed * Time.deltaTime);
        if (Input.GetKey(KeyCode.E)) Rotate(-rotateSpeed * Time.deltaTime);

        rocketOrbit.CalculatePositionVelocityatTime(Time.time);
        transform.position = rocketOrbit.position;

        rocketOrbit.DrawOrbitalLine(lineRenderer, orbitalLineNumberOfPoints);

        UIController.SetText(rocketOrbit.ToString());

        Debug.DrawLine(transform.position, transform.position+GetDirection());
        rocketOrbit.CalculateExtraVariables();
    }


    void Rotate(float angleRad) {
        transform.Rotate(new Vector3(0, 0, Mathf.Rad2Deg*angleRad));
    }

    Vector3 GetDirection() {
        float dirangle = Mathf.Deg2Rad*(transform.rotation.eulerAngles.z+90);
        //Debug.Log(dirangle);
        float x = Mathf.Cos(dirangle);
        float y = Mathf.Sin(dirangle);

        return new Vector3(x, y);
    }

}
