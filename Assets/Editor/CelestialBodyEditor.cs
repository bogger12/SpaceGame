using UnityEngine;
using UnityEditor;


[CustomEditor(typeof(CelestialBody))]
[CanEditMultipleObjects]
public class CelestialBodyEditor : Editor
{


    public bool rotateFoldout = true;
    public bool orbitFoldout = true;
    public bool lineFoldout = true;

    SerializedProperty mass;
    SerializedProperty Radius;
    SerializedProperty hasCustomSOI;
    SerializedProperty radiusSOI;

    SerializedProperty rotate;
    SerializedProperty rotatePeriod;

    SerializedProperty hasOrbit;
    SerializedProperty customOrbit;

    SerializedProperty e;
    SerializedProperty a;
    SerializedProperty w;
    SerializedProperty i;
    SerializedProperty omega;
    SerializedProperty f;

    SerializedProperty bodyOfInfluence;
    SerializedProperty initVelocity;

    SerializedProperty hasOrbitalLine;
    SerializedProperty lineColor;
    SerializedProperty lineWidth;


    private void OnEnable()
    {
        mass = serializedObject.FindProperty("mass");
        Radius = serializedObject.FindProperty("radius");
        hasCustomSOI = serializedObject.FindProperty("hasCustomSOI");
        radiusSOI = serializedObject.FindProperty("radiusSOI");

        rotate = serializedObject.FindProperty("rotate");
        rotatePeriod = serializedObject.FindProperty("rotatePeriod");

        hasOrbit = serializedObject.FindProperty("hasOrbit");
        customOrbit = serializedObject.FindProperty("customOrbit");

        e = serializedObject.FindProperty("e");
        a = serializedObject.FindProperty("a");
        w = serializedObject.FindProperty("w");
        i = serializedObject.FindProperty("i");
        omega = serializedObject.FindProperty("omega");
        f = serializedObject.FindProperty("f");

        bodyOfInfluence = serializedObject.FindProperty("bodyOfInfluence");
        initVelocity = serializedObject.FindProperty("initVelocity");

        hasOrbitalLine = serializedObject.FindProperty("hasOrbitalLine");
        lineColor = serializedObject.FindProperty("lineColor");
        lineWidth = serializedObject.FindProperty("lineWidth");
    }


    public override void OnInspectorGUI() {
        // If we call base the default inspector will get drawn too.
        // Remove this line if you don't want that to happen.

        //base.OnInspectorGUI();

        // This gets the current values from all serialized fields into the serialized "clone"
        serializedObject.Update();


        CelestialBody body = target as CelestialBody;

        Undo.RecordObject(body, "Change Celestial Body fields");


        mass.floatValue = EditorGUILayout.FloatField("Mass:", mass.floatValue);
        Radius.floatValue = EditorGUILayout.FloatField("Radius:", Radius.floatValue);
        body.Radius = Radius.floatValue;
        EditorGUILayout.BeginHorizontal();
        hasCustomSOI.boolValue = EditorGUILayout.Toggle("Override SOI", hasCustomSOI.boolValue);
        GUI.enabled = hasCustomSOI.boolValue;
        radiusSOI.floatValue = EditorGUILayout.FloatField(radiusSOI.floatValue);
        GUI.enabled = true;
        EditorGUILayout.EndHorizontal();


        EditorGUILayout.BeginHorizontal();
        rotate.boolValue = EditorGUILayout.Toggle("Rotation Period (s)", rotate.boolValue);
        GUI.enabled = rotate.boolValue;
        rotatePeriod.floatValue = EditorGUILayout.FloatField(rotatePeriod.floatValue);
        GUI.enabled = true;
        EditorGUILayout.EndHorizontal();

        EditorGUILayout.BeginHorizontal();
        orbitFoldout = EditorGUILayout.Foldout(orbitFoldout, "Orbit", true);
        hasOrbit.boolValue = EditorGUILayout.Toggle(hasOrbit.boolValue);
        EditorGUILayout.EndHorizontal();

        if (orbitFoldout) {
            GUI.enabled = hasOrbit.boolValue;
            customOrbit.boolValue = EditorGUILayout.Toggle("Custom Orbit", customOrbit.boolValue);
            if (customOrbit.boolValue) {
                float initwidth = EditorGUIUtility.labelWidth;
                EditorGUIUtility.labelWidth = 220;
                e.floatValue = EditorGUILayout.FloatField("Eccentricity", e.floatValue);
                a.floatValue = EditorGUILayout.FloatField("Semi-Major Axis (m)", a.floatValue);
                w.floatValue = EditorGUILayout.FloatField("Argument of Periapsis (rad)", w.floatValue);
                i.floatValue = EditorGUILayout.FloatField("Inclination (rad)", i.floatValue);
                omega.floatValue = EditorGUILayout.FloatField("Longtitude of Right-Ascending Node (rad)", omega.floatValue);
                f.floatValue = EditorGUILayout.FloatField("True Anomaly (rad)", f.floatValue);
                EditorGUIUtility.labelWidth = initwidth;
                EditorGUILayout.LabelField("Body of Influence:");
                bodyOfInfluence.objectReferenceValue = (CelestialBody)EditorGUILayout.ObjectField(bodyOfInfluence.objectReferenceValue, typeof(CelestialBody), true);
            }
            else {
                initVelocity.vector3Value = EditorGUILayout.Vector3Field("Initial Velocity:", initVelocity.vector3Value);
                EditorGUILayout.LabelField("Body of Influence:");
                bodyOfInfluence.objectReferenceValue = (CelestialBody)EditorGUILayout.ObjectField(bodyOfInfluence.objectReferenceValue, typeof(CelestialBody), true);
            }
            GUI.enabled = true;
        }

        EditorGUILayout.BeginHorizontal();
        lineFoldout = EditorGUILayout.Foldout(lineFoldout, "Orbital Line", true);
        hasOrbitalLine.boolValue = EditorGUILayout.Toggle(hasOrbitalLine.boolValue);
        EditorGUILayout.EndHorizontal();

        if (lineFoldout) {
            GUI.enabled = hasOrbitalLine.boolValue;
            lineColor.colorValue = EditorGUILayout.ColorField("Line Color:", lineColor.colorValue);
            lineWidth.floatValue = EditorGUILayout.FloatField("Line Width:", lineWidth.floatValue);
            GUI.enabled = true;
        }

        serializedObject.ApplyModifiedProperties();
        //save and load buttons
        if (EditorApplication.isPlaying) {
            if (GUILayout.Button("Apply")) {
                body.Start();
            }
        }
        Debug.DrawCircle(body.transform.position, body.CalculateSphereOfInfluence(), 20, Color.red, 100f);

    }
}