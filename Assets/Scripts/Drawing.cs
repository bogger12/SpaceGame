using UnityEngine;
using System.Collections;

public class Drawing
{
    public static void DrawCircle(Vector3 position, float radius, int segments, Color color, float duration) {
        GameObject surrogate = new GameObject("DEBUG");
        LineRenderer lines = surrogate.AddComponent<LineRenderer>();
        lines.positionCount = segments;
        lines.startColor = color;
        lines.endColor = color;
        lines.material = new Material(Shader.Find("Sprites/Default"));
        lines.sortingLayerID = SortingLayer.NameToID("Background");
        lines.loop = true;
        lines.startWidth = 0.1f;
        lines.endWidth = 0.1f;

        // If either radius or number of segments are less or equal to 0, skip drawing
        if (radius <= 0.0f || segments <= 0) {
            return;
        }

        // Single segment of the circle covers (360 / number of segments) degrees
        float angleStep = (2*Mathf.PI / segments);

        // lineStart and lineEnd variables are declared outside of the following for loop
        Vector3 point = Vector3.zero;

        for (int i = 0; i < segments; i++) {
            // Line start is defined as starting angle of the current segment (i)
            point.x = Mathf.Cos(angleStep * i);
            point.y = Mathf.Sin(angleStep * i);

            // Results are multiplied so they match the desired radius
            point *= radius;

            // Results are offset by the desired position/origin 
            point += position;

            // Points are connected using DrawLine method and using the passed color
            //Debug.DrawLine(lineStart, lineEnd, color, duration);
            lines.SetPosition(i, point);
        }
        GameObject.Destroy(surrogate, duration);
    }
}

