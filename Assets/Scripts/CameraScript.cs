using System.Collections;
using System.Collections.Generic;
using UnityEngine;



public class CameraScript : MonoBehaviour
{

    public Parallax parallaxObject;

    public GameObject rocket;

    [Range(0,10)]
    public float moveSpeed;


    private Vector3 roughPos;

    void Start()
    {
        roughPos = transform.position;
    }

    // Update is called once per frame
    void Update()
    {
        if (Input.GetKey(KeyCode.I)) Move(Vector3.up * moveSpeed * Time.deltaTime);
        if (Input.GetKey(KeyCode.K)) Move(Vector3.down * moveSpeed * Time.deltaTime);
        if (Input.GetKey(KeyCode.J)) Move(Vector3.left * moveSpeed * Time.deltaTime);
        if (Input.GetKey(KeyCode.L)) Move(Vector3.right * moveSpeed * Time.deltaTime);

        //transform.position = VPixelSnap(roughPos);
        SetPosition(VPixelSnap(rocket.transform.position));

        parallaxObject.SetParallaxPoint(transform.position);
    }

    void Move(Vector3 vmove) {
        roughPos += vmove;
    }

    Vector3 VPixelSnap(Vector3 v) {
        float pixelSnap(float a) { return GameSystem.pixelUnit * Mathf.Round(a / GameSystem.pixelUnit); }
        return new Vector3(pixelSnap(v.x), pixelSnap(v.y), pixelSnap(v.z));
    }

    void SetPosition(Vector2 v) {
        float initZ = transform.position.z;
        transform.position = v;
        transform.position += new Vector3(0, 0, initZ);

    }

}
