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

    [Header("Input Values:")]
    [Range(0, 10)]
    public float moveSpeed;
    [Range(0f, 1f)]
    public float zoomSpeed;

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

    // Update is called once per frame
    void Update()
    {
        if (Input.GetKey(KeyCode.I)) cameraScript.Move(moveSpeed * Time.deltaTime * Vector3.up);
        if (Input.GetKey(KeyCode.K)) cameraScript.Move(moveSpeed * Time.deltaTime * Vector3.down);
        if (Input.GetKey(KeyCode.J)) cameraScript.Move(moveSpeed * Time.deltaTime * Vector3.left);
        if (Input.GetKey(KeyCode.L)) cameraScript.Move(moveSpeed * Time.deltaTime * Vector3.right);

        if (Input.GetMouseButtonDown(1)) {
            cameraScript.setCenterOn(null);
            initMouse = screenCameraComponent.ScreenToWorldPoint(Input.mousePosition);
            lastMouse = initMouse;
            initCameraPos = cameraScript.roughPos;
        }
        else if (Input.GetMouseButton(1)) {
            Vector2 worlddistance = screenCameraComponent.ScreenToWorldPoint(Input.mousePosition) - lastMouse;
            lastMouse = screenCameraComponent.ScreenToWorldPoint(Input.mousePosition);
            cameraScript.Move(-worlddistance * GameSystem.screenScale);
        }
        else if (Input.GetMouseButtonUp(1)) {
            Vector2 worlddistance = screenCameraComponent.ScreenToWorldPoint(Input.mousePosition) - initMouse;
            cameraScript.roughPos = initCameraPos;
            cameraScript.Move(-worlddistance * GameSystem.screenScale);
        }

        // Use scrollwheel to change scale
        cameraScript.Scale(1f - zoomSpeed * Input.mouseScrollDelta.y);

        // Set the cursor position to the relative mouse position
        cursor.transform.localPosition = GameSystem.VPixelSnap((GameSystem.V3SetZ((screenCameraComponent.ScreenToWorldPoint(Input.mousePosition)), cursor.transform.localPosition.z)));
        // No need to scale the pixels to snap to when changing scale, as it is a parent of the camera gameObject.

        if (Input.GetMouseButtonDown(0)) {
            Vector3 clickPos = screenCameraComponent.ScreenToWorldPoint(Input.mousePosition);
            clickPos = GameSystem.V3SetZ((clickPos * GameSystem.screenScale + cameraScript.roughPos), 0.5f);

            RaycastHit2D raycasthit = Physics2D.Raycast(clickPos, Vector3.forward);
            if (raycasthit.collider!=null) { 
                if (raycasthit.collider.CompareTag(GameSystem.celestialBodyTag) || raycasthit.collider.CompareTag(GameSystem.rocketTag)) {
                    cameraScript.setCenterOn(raycasthit.transform);
                }
            }
        }

    }
}

