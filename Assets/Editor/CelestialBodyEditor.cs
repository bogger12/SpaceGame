using UnityEngine;
using UnityEditor;


[CustomEditor(typeof(CelestialBody))]
[CanEditMultipleObjects]
public class CelestialBodyEditor : Editor
{


    public bool rotateFoldout = true;
    public bool orbitFoldout = true;
    public bool lineFoldout = true;


    public override void OnInspectorGUI() {
        // If we call base the default inspector will get drawn too.
        // Remove this line if you don't want that to happen.

        //base.OnInspectorGUI();

        CelestialBody body = target as CelestialBody;

        body.mass = EditorGUILayout.FloatField("Mass:", body.mass);
        body.Radius = EditorGUILayout.Slider("Radius:", body.Radius, 0f, 50f);

        EditorGUILayout.BeginHorizontal();
        rotateFoldout = EditorGUILayout.Foldout(rotateFoldout, "Rotation", true);
        body.rotate = EditorGUILayout.Toggle(body.rotate);
        EditorGUILayout.EndHorizontal();


        if (rotateFoldout) {
            GUI.enabled = body.rotate;
            body.rotateSpeed = EditorGUILayout.Slider(body.rotateSpeed, 0f, 10f);
            GUI.enabled = true;
        }

        EditorGUILayout.BeginHorizontal();
        orbitFoldout = EditorGUILayout.Foldout(orbitFoldout, "Orbit", true);
        body.hasOrbit = EditorGUILayout.Toggle(body.hasOrbit);
        EditorGUILayout.EndHorizontal();

        if (orbitFoldout) {
            GUI.enabled = body.hasOrbit;
            body.customOrbit = EditorGUILayout.Toggle("Custom Orbit", body.customOrbit);
            if (body.customOrbit) {
                float initwidth = EditorGUIUtility.labelWidth;
                EditorGUIUtility.labelWidth = 220;
                body.e = EditorGUILayout.FloatField("Eccentricity", body.e);
                body.a = EditorGUILayout.FloatField("Semi-Major Axis (m)", body.a);
                body.w = EditorGUILayout.FloatField("Argument of Periapsis (rad)", body.w);
                body.i = EditorGUILayout.FloatField("Inclination (rad)", body.i);
                body.omega = EditorGUILayout.FloatField("Longtitude of Right-Ascending Node (rad)", body.omega);
                body.f = EditorGUILayout.FloatField("True Anomaly (rad)", body.f);
                EditorGUIUtility.labelWidth = initwidth;
                EditorGUILayout.LabelField("Body of Influence:");
                body.bodyOfInfluence = (GameObject)EditorGUILayout.ObjectField(body.bodyOfInfluence, typeof(GameObject), true);
                //save and load buttons
                if (EditorApplication.isPlaying) { 
                    if (GUILayout.Button("Apply")) {
                        body.Start();
                    }
                }
            }
            else {
                body.initVelocity = EditorGUILayout.Vector3Field("Initial Velocity:", body.initVelocity);
                EditorGUILayout.LabelField("Body of Influence:");
                body.bodyOfInfluence = (GameObject)EditorGUILayout.ObjectField(body.bodyOfInfluence, typeof(GameObject), true);
            }
            GUI.enabled = true;
        }

        EditorGUILayout.BeginHorizontal();
        lineFoldout = EditorGUILayout.Foldout(lineFoldout, "Orbital Line", true);
        body.hasOrbitalLine = EditorGUILayout.Toggle(body.hasOrbitalLine);
        EditorGUILayout.EndHorizontal();

        if (lineFoldout) {
            GUI.enabled = body.hasOrbitalLine;
            body.lineColor = EditorGUILayout.ColorField("Line Color:", body.lineColor);
            body.lineWidth = EditorGUILayout.FloatField("Line Width:", body.lineWidth);
            GUI.enabled = true;
        }
    }
}