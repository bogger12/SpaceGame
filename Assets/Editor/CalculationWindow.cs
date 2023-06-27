using UnityEditor;
using UnityEngine;
using UnityEngine.UIElements;
using UnityEditor.UIElements;


public class CalculationWindow : EditorWindow {
    [MenuItem("Window/UI Toolkit/CalculationWindow")]
    public static void ShowExample() {
        CalculationWindow wnd = GetWindow<CalculationWindow>();
        wnd.titleContent = new GUIContent("CalculationWindow");
    }

    public void CreateGUI() {
        // Each editor window contains a root VisualElement object
        VisualElement root = rootVisualElement;

        // VisualElements objects can contain other VisualElement following a tree hierarchy.


        // find velocity needed, given gravity, masses, and radius

        root.Add(new Label("\nInput values to find velocity needed:"));

        
        FloatField radius = new FloatField("Radius");
        root.Add(radius);

        ObjectField celestialBodyField = new ObjectField("Celestial Body");
        celestialBodyField.objectType = typeof(CelestialBody);
        root.Add(celestialBodyField);

        Button but = new Button();
        but.text = "GO";
        Label resultLabel = new Label("output is here");
        resultLabel.text = "";


        but.RegisterCallback<ClickEvent>((ClickEvent evt) => {
            if (celestialBodyField.value is null) resultLabel.text = "Make sure all fields are filled";
            CelestialBody bodyOfInfluence = (CelestialBody) celestialBodyField.value;
            // To find centripetal velocity:
            float velocity = Mathf.Sqrt((Orbit.gravconst * bodyOfInfluence.mass)/(radius.value));
            resultLabel.text += "Velocity required: " + velocity + "m/s";
        });
        root.Add(but);
        root.Add(resultLabel);

    }
}