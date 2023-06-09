using System.Collections;
using System.Collections.Generic;
using UnityEngine;



public class CameraScript : MonoBehaviour
{

    public Parallax parallaxObject;

    public Transform centerOn = null;

    public Vector3 roughPos;

    private Camera camera;

    void Start()
    {
        roughPos = transform.position;
        camera = GetComponent<Camera>();
    }

    // LateUpdate is used as we want to track the object after it has moved (changed order in Script Execution Order so no longer needed)
    void Update()
    {
        if (centerOn != null) { roughPos = GameSystem.V3SetZ(centerOn.position, roughPos.z); }
        //transform.position = GameSystem.VPixelSnap(roughPos);
        transform.position = roughPos;

        parallaxObject.SetParallaxPoint(transform.position);
    }

    public void Move(Vector3 vmove) {
        roughPos += vmove;
    }

    public void Scale(float multiply) {
        Vector3 multiplyOnlyXY(Vector3 v, float multiply) {
            return (Vector3)((Vector2)v * multiply) + Vector3.forward;
        }
        camera.orthographicSize *= multiply;
        transform.localScale = multiplyOnlyXY(transform.localScale, multiply);
        GameSystem.screenScale *= multiply;
    }

    public void setCenterOn(Transform newcenter) {
        if (newcenter==null && centerOn != null) {
            roughPos = GameSystem.V3SetZ(centerOn.position, roughPos.z);
        }
        centerOn = newcenter;
    }

}
