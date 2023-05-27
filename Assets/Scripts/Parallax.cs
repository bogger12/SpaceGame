using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Parallax : MonoBehaviour
{

    public GameObject[] parallaxObjects;
    public float[] parallaxValues;

    public float objectScale = 1f;

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
            g.transform.localPosition = (parallaxPoint) - pvalue*(objectScale)*(parallaxPoint);
            g.transform.localScale = Vector2.one * objectScale;
            //Debug.Log(parallaxPoint * objectScale);
        }
    }

}
