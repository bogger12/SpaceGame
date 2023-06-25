using UnityEngine;
using System.Collections;
using UnityEngine.LowLevel;
using System.Text;

class CustomLoop : MonoBehaviour {

    [RuntimeInitializeOnLoadMethod]
    public static void AppStart() {
        if (!Application.isPlaying) return;
        PlayerLoopSystem defaultSystems = PlayerLoop.GetDefaultPlayerLoop();

        var OrbitUpdate = new PlayerLoopSystem() {
            updateDelegate = OrbitUpdates,
            type = typeof(CustomLoop)
        };
        InsertSubSystem<UnityEngine.PlayerLoop.Update.ScriptRunBehaviourUpdate>(ref defaultSystems, OrbitUpdate);
        PlayerLoop.SetPlayerLoop(defaultSystems);

        //var sb = new StringBuilder();
        //RecursivePlayerLoopPrint(defaultSystems, sb, 0);
        //Debug.Log(sb.ToString());
    }


    // Inserts a custom subsystem at the start of the array that system is in
    private static bool InsertSubSystem<T>(ref PlayerLoopSystem system, PlayerLoopSystem replacement) {
        if (system.type == typeof(T)) {
            //system = replacement;
            return true;
        }
        if (system.subSystemList != null) {
            for (var i = 0; i < system.subSystemList.Length; i++) {
                if (InsertSubSystem<T>(ref system.subSystemList[i], replacement)) {
                    PlayerLoopSystem[] newLoopSystem = new PlayerLoopSystem[system.subSystemList.Length + 1];
                    newLoopSystem[0] = replacement;
                    for (int k=1;k < system.subSystemList.Length+1;k++) {
                        newLoopSystem[k] = system.subSystemList[k - 1];
                    }
                    system.subSystemList = newLoopSystem;
                    return false;
                }
            }
        }
        return false;
    }

    //private static void RecursivePlayerLoopPrint(PlayerLoopSystem def, StringBuilder sb, int depth) {
    //    if (depth == 0) {
    //        sb.AppendLine("ROOT NODE");
    //    }
    //    else if (def.type != null) {
    //        for (int i = 0; i < depth; i++) {
    //            sb.Append("\t");
    //        }
    //        sb.AppendLine(def.type.Name);
    //    }
    //    if (def.subSystemList != null) {
    //        depth++;
    //        foreach (var s in def.subSystemList) {
    //            RecursivePlayerLoopPrint(s, sb, depth);
    //        }
    //        depth--;
    //    }
    //}

    public static void OrbitUpdates() {
        GameObject.FindObjectOfType<RocketScript>().OrbitalUpdate();
        Transform rootBody = GameObject.FindGameObjectWithTag("Root Body").transform;
        RecursiveOrbitalUpdate(rootBody); // Calls each OrbitalUpdate function in recursive order from the root
    }

    public static void RecursiveOrbitalUpdate(Transform rootBody) {
        rootBody.GetComponent<CelestialBody>().OrbitalUpdate();
        if (rootBody.transform.childCount == 0) return;
        else {
            foreach(Transform t in rootBody) {
                RecursiveOrbitalUpdate(t);
            }
        }
    }
}

