using UnityEngine;
using System.Collections;
using Unity.VisualScripting;
using static UnityEngine.UI.CanvasScaler;

public class InputManager : MonoBehaviour
{
    public GameObject rocket;
    public GameObject camera;
    public GameObject screenCamera;
    public GameObject cursor;

    private Texture2D cursorTex;

    private RocketScript rocketScript;
    private CameraScript cameraScript;
    private Camera screenCameraComponent;

    private Vector3 initCameraPos;
    private Vector3 initMouse;
    private Vector3 lastMouse;

    // Use this for initialization
    void Start()
    {
        rocketScript = rocket.GetComponent<RocketScript>();
        cameraScript = camera.GetComponent<CameraScript>();
        screenCameraComponent = screenCamera.GetComponent<Camera>();
        cursorTex = cursor.GetComponent<SpriteRenderer>().sprite.texture;

        Cursor.visible = false;
    }

    Vector3 vXY(Vector3 v, float z) {
        Vector3 result = (Vector3)((Vector2)v);
        return result + new Vector3(0, 0, z);
    }

    // Update is called once per frame
    void Update()
    {
        if (Input.GetKey(KeyCode.I)) cameraScript.Move(cameraScript.moveSpeed * Time.deltaTime * Vector3.up);
        if (Input.GetKey(KeyCode.K)) cameraScript.Move(cameraScript.moveSpeed * Time.deltaTime * Vector3.down);
        if (Input.GetKey(KeyCode.J)) cameraScript.Move(cameraScript.moveSpeed * Time.deltaTime * Vector3.left);
        if (Input.GetKey(KeyCode.L)) cameraScript.Move(cameraScript.moveSpeed * Time.deltaTime * Vector3.right);

        if (Input.GetKey(KeyCode.O)) cameraScript.Scale(1f - cameraScript.zoomSpeed * Time.deltaTime);
        if (Input.GetKey(KeyCode.P)) cameraScript.Scale(1f + cameraScript.zoomSpeed * Time.deltaTime);

        if (Input.GetMouseButtonDown(0)) {
            initMouse = screenCameraComponent.ScreenToWorldPoint(Input.mousePosition);
            lastMouse = initMouse;
            initCameraPos = cameraScript.roughPos;
        }
        else if (Input.GetMouseButton(0)) {
            Vector2 worlddistance = screenCameraComponent.ScreenToWorldPoint(Input.mousePosition) - lastMouse;
            lastMouse = screenCameraComponent.ScreenToWorldPoint(Input.mousePosition);
            cameraScript.Move(-worlddistance * cameraScript.scale);
        }
        else if (Input.GetMouseButtonUp(0)) {
            Vector2 worlddistance = screenCameraComponent.ScreenToWorldPoint(Input.mousePosition) - initMouse;
            cameraScript.roughPos = initCameraPos;
            cameraScript.Move(-worlddistance*cameraScript.scale);
        }

        //cursor.transform.position = vXY((screenCameraComponent.ScreenToWorldPoint(Input.mousePosition) * cameraScript.scale + cameraScript.roughPos), cursor.transform.position.z);
        //cursor.transform.position += (new Vector3(+cursorTex.width * GameSystem.pixelUnit / 2, -cursorTex.height * GameSystem.pixelUnit / 2));

        cursor.transform.localPosition = GameSystem.VSnap((vXY((screenCameraComponent.ScreenToWorldPoint(Input.mousePosition)), cursor.transform.localPosition.z)), GameSystem.pixelUnit);
    }
}

