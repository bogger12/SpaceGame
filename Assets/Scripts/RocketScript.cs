using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RocketScript : MonoBehaviour {

    public UIController UIController;
    public LineRenderer lineRenderer;
    public GameObject bodyOfInfluence;

    public Orbit rocketOrbit;

    private float mass = 1f;



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

        if (Input.GetKey(KeyCode.Q)) Rotate(1f * Time.deltaTime);
        if (Input.GetKey(KeyCode.E)) Rotate(-1f * Time.deltaTime);

        rocketOrbit.CalculatePositionVelocityatTime(Time.time);
        transform.position = rocketOrbit.position;

        rocketOrbit.DrawOrbitalLine(lineRenderer, 50);

        UIController.SetText(rocketOrbit.ToString());

        Debug.DrawLine(transform.position, transform.position+GetDirection());
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
