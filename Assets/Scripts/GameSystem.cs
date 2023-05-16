using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class GameSystem
{
    public static int PPU = 16;
    public static float pixelUnit { get { return 1 / (float)PPU; } }

    public static bool DEBUG = true;


    public static void Rotate(Transform transform, float angleRad) {
        transform.Rotate(new Vector3(0, 0, Mathf.Rad2Deg * angleRad));
    }
}
