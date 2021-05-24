using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System;
using System.Globalization;

public class Coil:MonoBehaviour
{
    private static CultureInfo Culture = CultureInfo.InvariantCulture;
    
    public int windingNum = 2; // число витков
    public int turnsNum = 2; // число намоток катушки
    public float wireWidth = 0.1f; // м, ширина намотки
    public double i = 1.0f; // А, величина силы тока, протекающего по контуру
    public double R = 1.0f; // м, внешний радиус контура S
    public GameObject ringPrefab;
    
    public InputField ui_windingNum;
    public InputField ui_turnsNum;
    public InputField ui_wireWidth;
    public InputField ui_i;
    public InputField ui_R; 
 //   public GameObject ui_ringPrefab;
    
    [HideInInspector]
    public List<Transform> rings = new List<Transform>();
    
    void Awake()
    {
        StateController.ESetState += OnSetState;
    }

    void OnDestroy()
    {
        StateController.ESetState -= OnSetState;
    }

	public void OnMouseUpAsButton()
	{
		Debug.Log("Clicked");
	}

    float ConvertFromSIToUnitySize(float value)
    {
        float newValue = value * 10f;
        return newValue;
    }
    
    float ConvertFromUnitySizeToSI(float value)
    {
        float newValue = value / 10f;
        return newValue;
    }

    public void SpawnRings()
    {
        float halfLength = windingNum / 2f * (10 * wireWidth) - (10 * wireWidth)/2;
        //float halfLength = windingNum * wireWidth / 2f;
        
        // Переходим к размерности Unity
        Vector3 startPos = transform.position - transform.forward * halfLength;

        for (int x = 0; x < windingNum; x++)
        {
            var clone = Instantiate(ringPrefab, startPos + transform.forward * (10 * wireWidth) * x,
                transform.rotation);
            clone.transform.localScale = new Vector3((float) (10f * R), (float) (10f * R), (10f * wireWidth));
            rings.Add(clone.transform);
            clone.transform.SetParent(this.transform, true);
        }
    }

    public void Reinit()
    {
        if (gameObject.activeInHierarchy)
            SpawnRings();
        
        gameObject.GetComponent<MeshRenderer>().enabled = false;
    }

    void OnSetState(SceneState state)
    {
        if (state==SceneState.Setup)
        {
            RemoveRings();
        }
    }

    void RemoveRings()
    {
        foreach (var r in rings)
            Destroy(r.gameObject);

        rings.Clear();
        
        gameObject.GetComponent<MeshRenderer>().enabled = true;
    }

    public void SetUIValues()
    {
        ui_windingNum.text = windingNum.ToString();
        ui_turnsNum.text = turnsNum.ToString();
        ui_wireWidth.text = wireWidth.ToString();
        ui_i.text = i.ToString();
        ui_R.text = R.ToString();
    }

    public void GetUIValues()
    {
        if (!(ui_windingNum is null))windingNum = Convert.ToInt32(ui_windingNum.text, Culture);
        if (!(ui_turnsNum is null))turnsNum = Convert.ToInt32(ui_turnsNum.text, Culture);
        if (!(ui_wireWidth is null))wireWidth = Convert.ToSingle(ui_wireWidth.text, Culture);
        
        Double doubleNum;
        if (!(ui_i is null) && Double.TryParse(ui_i.text, out doubleNum)) i = doubleNum;
        if (!(ui_R is null) && Double.TryParse(ui_R.text, out doubleNum)) R = doubleNum;
    }
}