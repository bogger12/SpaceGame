using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class GameSystem
{
    public static int screenWidth = 512;
    public static int screenHeight = 448;

    public static int PPU = 16;
    public static float pixelUnit { get { return 1 / (float)PPU; } }

    public static bool DEBUG = true;
    public static bool LOG_ELEMENTS = false;

    public static string celestialBodyTag = "Celestial Body";


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
}
