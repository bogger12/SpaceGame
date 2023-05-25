using System.Collections;
using System.Collections.Generic;
using UnityEngine;



public class CameraScript : MonoBehaviour
{

    public Parallax parallaxObject;
    [Range(0,10)]
    public float moveSpeed;
    [Range(0, 10)]
    public float zoomSpeed;

    public Transform centerOn = null;

    public Vector3 roughPos;

    private Camera camera;

    public float scale = 1f;

    void Start()
    {
        roughPos = transform.position;
        camera = GetComponent<Camera>();
    }

    // LateUpdate is used as we want to track the object after it has moved
    void LateUpdate()
    {

        if (centerOn != null) {
            SetPosition(GameSystem.VPixelSnap(centerOn.position));
        }
        else {
            transform.position = GameSystem.VPixelSnap(roughPos);
        }

        parallaxObject.SetParallaxPoint(transform.position);
    }

    public void Move(Vector3 vmove) {
        roughPos += vmove;
    }

    public void SetPosition(Vector2 v) {
        float initZ = transform.position.z;
        transform.position = new Vector3(v.x, v.y, initZ);
    }

    public void Scale(float multiply) {
        Vector3 multiplyOnlyXY(Vector3 v, float multiply) {
            return (Vector3)((Vector2)v * multiply) + Vector3.forward;
        }
        camera.orthographicSize *= multiply;
        transform.localScale = multiplyOnlyXY(transform.localScale, multiply);
        scale *= multiply;
        parallaxObject.objectScale = scale;
    }

}
