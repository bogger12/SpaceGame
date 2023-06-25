using UnityEditor;
using UnityEngine;
using UnityEngine.UIElements;
using UnityEditor.UIElements;


public class TestWindow : EditorWindow
{
    [MenuItem("Window/UI Toolkit/TestWindow")]
    public static void ShowExample()
    {
        TestWindow wnd = GetWindow<TestWindow>();
        wnd.titleContent = new GUIContent("TestWindow");
    }

    public void CreateGUI()
    {
        // Each editor window contains a root VisualElement object
        VisualElement root = rootVisualElement;

        // VisualElements objects can contain other VisualElement following a tree hierarchy.

        root.Add(new Label("\nInput values for orbit testing"));
        Vector3Field inputR = new Vector3Field("Position");
        Vector3Field inputV = new Vector3Field("Velocity");
        ObjectField celestialBodyField = new ObjectField("Celestial Body");
        celestialBodyField.objectType = typeof(CelestialBody);
        root.Add(inputR);
        root.Add(inputV);
        root.Add(celestialBodyField);

        Button but = new Button();
        Label resultLabel = new Label("output is here");
        but.text = "GO";
        but.RegisterCallback<ClickEvent>((ClickEvent evt) => {
            if (celestialBodyField.value is null) resultLabel.text = "Make sure all fields are filled";
            Orbit testOrbit = new Orbit(inputR.value, inputV.value, (CelestialBody)celestialBodyField.value, 1);
            testOrbit.CalculatePositionVelocityatTime(Time.time, true, false);
            resultLabel.text = testOrbit.ToString();
        });
        root.Add(but);
        root.Add(resultLabel);

    }
}