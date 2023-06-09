using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Parallax : MonoBehaviour
{

    public GameObject[] parallaxObjects;
    public float[] parallaxValues;

    //[HideInInspector]
    public Vector2 parallaxPoint;

    // Update is called once per frame
    void Update()
    {
        
    }

    public void SetParallaxPoint(Vector2 point) {
        parallaxPoint = point;
        SetTransformsBasedOnPValues();
    }


    void SetTransformsBasedOnPValues() {
        GameObject g;
        float pvalue;
        for (int i = 0; i < parallaxObjects.Length; i++) {
            g = parallaxObjects[i];
            pvalue = parallaxValues[i];
            //g.transform.localPosition = pvalue * (parallaxPoint);
            //g.transform.localPosition = GameSystem.VPixelSnap((parallaxPoint) - pvalue*(GameSystem.screenScale) *(parallaxPoint));
            g.transform.localPosition = (parallaxPoint) - pvalue*(GameSystem.screenScale) *(parallaxPoint);
            g.transform.localScale = Vector2.one * GameSystem.screenScale;
            //Debug.Log(parallaxPoint * objectScale);
        }
    }

}
