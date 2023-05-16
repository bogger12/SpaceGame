using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class GameSystem
{
    public static int PPU = 16;
    public static float pixelUnit {
        get { return 1 / (float)PPU; }
    }


}
