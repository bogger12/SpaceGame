using System.IO;
using UnityEditor;
using UnityEngine;
using UnityEngine.UIElements;
using UnityEditor.UIElements;

public class CreateSystem : EditorWindow {

    string csvPath;

    [MenuItem("Window/UI Toolkit/CreateSystem")]
    public static void ShowExample() {
        CreateSystem wnd = GetWindow<CreateSystem>();
        wnd.titleContent = new GUIContent("CreateSystem");
    }

    public void CreateGUI() {
        // Each editor window contains a root VisualElement object
        VisualElement root = rootVisualElement;

        Undo.RecordObject(this, "Change CreateSystem fields");

        // VisualElements objects can contain other VisualElement following a tree hierarchy.

        root.Add(new Label("\nInput CSV file for planets"));

        TextField csvPathField = new TextField("Enter file path:");
        csvPathField.value = "Assets/Resources/";
        root.Add(csvPathField);

        TextField spritePathField = new TextField("Enter Path to Default Sprite:");
        //spritePathField.value = "";
        root.Add(spritePathField);


        FloatField massMulti = new FloatField("Mass Multiplier");
        massMulti.value = 1E-6f;
        root.Add(massMulti);

        FloatField distanceMulti = new FloatField("Distance Multiplier");
        distanceMulti.value = 1E-6f;
        root.Add(distanceMulti);

        FloatField sizeMulti = new FloatField("Size Multiplier");
        sizeMulti.value = 1f;
        root.Add(sizeMulti);

        Button gobutton = new Button();
        gobutton.text = "Create System";
        root.Add(gobutton);

        Label resultLabel = new Label("output is here");
        root.Add(resultLabel);

        gobutton.RegisterCallback<ClickEvent>((ClickEvent evt) => {
            csvPath = csvPathField.value;
            if (csvPath is null) resultLabel.text = "Please enter path for CSV file";
            resultLabel.text = CreateSystemObjects(ReadCSV(csvPath, 21, 11), massMulti.value, distanceMulti.value, sizeMulti.value, spritePathField.value);
        });
    }

    public string[][] ReadCSV(string csvPath, int rows, int columns) {


        string removeChar(char removeC, string removeFrom) {
            while (removeFrom.Contains(removeC)) {
                for (int i = 0; i < removeFrom.Length; i++) {
                    if (removeFrom[i] == removeC) {
                        removeFrom = removeFrom.Substring(0, i) + removeFrom.Substring(i + 1, removeFrom.Length - (i + 1));
                        i++;
                    }
                }
            }
            return removeFrom;
        }


        string[][] returnstring = new string[rows][];
        StreamReader streamReader = new StreamReader(csvPath);

        bool EOF = false;

        int lineCount = 0;

        while (!EOF) {
            string currentLine = streamReader.ReadLine();
            if (currentLine == null) {
                EOF = true;
                break;
            }

            currentLine = removeChar('*', currentLine);
            while(currentLine.Contains('"')) {
                for (int i = 0; i < currentLine.Length; i++) {
                    if (currentLine[i] == '"') {
                        int commapos = currentLine.Substring(i).IndexOf(",");
                        int nexticommapos = currentLine.Substring(i+1).IndexOf('"')+1;

                        string b4ic = currentLine.Substring(0, i);
                        string b4c = currentLine.Substring(i + 1, commapos - 1);
                        string aftc = currentLine.Substring(i + commapos + 1, nexticommapos - commapos - 1);
                        string aftic = currentLine.Substring(i + nexticommapos + 1, currentLine.Length - (i + nexticommapos + 1));

                        currentLine = b4ic + b4c + aftc + aftic;
                        break;
                    }
                }
            }



            returnstring[lineCount++] = currentLine.Split(",");
        }
        return returnstring;

    }


    public string CreateSystemObjects(string[][] dataTable, float massMultiplier, float distanceMultiplier, float sizeMultiplier, string spritePath) {
        GameObject system = new GameObject("System");

        GameObject sun = new GameObject("Sun", typeof(SpriteRenderer), typeof(CircleCollider2D), typeof(CelestialBody));
        sun.transform.SetParent(system.transform);
        CelestialBody sunBody = sun.GetComponent<CelestialBody>();
        sun.GetComponent<SpriteRenderer>().sprite = Resources.Load<Sprite>("Sprites/Sun512") as Sprite;
        sun.GetComponent<SpriteRenderer>().drawMode = SpriteDrawMode.Sliced;

        sunBody.mass = 1.9885E+30f * massMultiplier;
        sunBody.Radius = 695700f * distanceMultiplier * sizeMultiplier;
        sunBody.density = 1408f * massMultiplier / (distanceMultiplier * distanceMultiplier * distanceMultiplier);
        sunBody.rotatePeriod = 609.12f; // * 60f; // multiply by 60 to get time in seconds
        sunBody.hasOrbit = false;
        sunBody.hasCustomSOI = true;
        sunBody.radiusSOI = 60000E+6f * distanceMultiplier;
        sun.transform.position = Vector2.zero;

        for (int i = 1; i < dataTable[0].Length;i++) {
            Debug.Log(dataTable[0][i]);
            GameObject planet = new GameObject(
                dataTable[0][i][0] + dataTable[0][i].Substring(1).ToLower(),
                typeof(SpriteRenderer),
                typeof(CircleCollider2D),
                typeof(CelestialBody)
            );
            planet.transform.SetParent(sun.transform);

            CelestialBody planetBody = planet.GetComponent<CelestialBody>();
            planet.GetComponent<SpriteRenderer>().sprite = Resources.Load<Sprite>("Sprites/"+spritePath) as Sprite;
            planet.GetComponent<SpriteRenderer>().drawMode = SpriteDrawMode.Sliced;

            planetBody.bodyOfInfluence = sun.GetComponent<CelestialBody>();

            planetBody.mass = (float.Parse(dataTable[1][i]) * 10E24f) * massMultiplier;

            planetBody.Radius = (float.Parse(dataTable[2][i]) / 2f) * distanceMultiplier * sizeMultiplier;

            planetBody.density = float.Parse(dataTable[3][i]) * massMultiplier / (distanceMultiplier * distanceMultiplier * distanceMultiplier);

            planetBody.rotatePeriod = float.Parse(dataTable[6][i]); // * 60f; // multiply by 60 to get time in seconds

            planetBody.a = (float.Parse(dataTable[9][i]) * 10E6f) * distanceMultiplier;
            planet.transform.position = new Vector2(planetBody.a, 0);

            planetBody.e = float.Parse(dataTable[14][i]);

            planetBody.lineColor = new Color(1, 1, 1, 0.3f);
        }
        return "Success!";
    }
}

