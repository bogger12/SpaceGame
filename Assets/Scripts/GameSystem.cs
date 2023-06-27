using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using static TMPro.SpriteAssetUtilities.TexturePacker_JsonArray;

public static class GameSystem
{
    public static int screenWidth = 512;
    public static int screenHeight = 448;

    public static float timeElapsed = 0f;
    public static float timeScale = 1f;

    public static int PPU = 16;
    public static float pixelUnit { get { return 1 / (float)PPU; } }

    public static float screenScale = 1f;

    public static bool DEBUG = true;
    public static bool LOG_ELEMENTS = true;

    public static string celestialBodyTag = "Celestial Body";
    public static string rocketTag = "Rocket";


    public static float CurrentTime() { return timeElapsed; }
    public static float DeltaTime() { return Time.deltaTime * timeScale; }
    public static float FixedDeltaTime() { return Time.fixedDeltaTime * timeScale; }

    public static void Rotate(Transform transform, float angleRad) {
        transform.Rotate(new Vector3(0, 0, Mathf.Rad2Deg * angleRad));
    }

    public static Vector3 VPixelSnap(Vector3 v) {
        float pixelSnap(float a) { return pixelUnit * Mathf.Round(a / pixelUnit); }
        return new Vector3(pixelSnap(v.x), pixelSnap(v.y), pixelSnap(v.z));
    }
    public static Vector3 VSnap(Vector3 v, float unit) {
        float pixelSnap(float a) { return unit * Mathf.Round(a / unit); }
        return new Vector3(pixelSnap(v.x), pixelSnap(v.y), pixelSnap(v.z));
    }

    public static Vector3 V3SetZ(Vector3 v, float z) {
        Vector3 result = (Vector3)((Vector2)v);
        return result + new Vector3(0, 0, z);
    }

    public static void TestAngle(Vector3 point, float initAngle, float angle) {
        Debug.DrawLine(point, point+(Vector3)PolarToCartesian(5f, initAngle), Color.blue);
        Debug.DrawLine(point, point+(Vector3)PolarToCartesian(5f, angle), Color.red);
    }

    public static Vector2 PolarToCartesian(float hyp, float theta) {
        return new Vector2(hyp * Mathf.Cos(theta), hyp * Mathf.Sin(theta));
    }

    public static Vector2 RotateVector(Vector2 v, float angle) {
        float x = Mathf.Cos(angle) * v.x - Mathf.Sin(angle) * v.y;
        float y = Mathf.Sin(angle) * v.x + Mathf.Cos(angle) * v.y;
        return new Vector2(x, y);
    }

    public static Vector3 AngleToDirection(float angle) {
        float dirangle = Mathf.Deg2Rad * (angle + 90);
        return new Vector3(Mathf.Cos(dirangle), Mathf.Sin(dirangle)).normalized;
    }

    public static void SetLineWidth(LineRenderer lr, float lineWidth) {
        lr.startWidth = lineWidth;
        lr.endWidth = lineWidth;
    }
}
