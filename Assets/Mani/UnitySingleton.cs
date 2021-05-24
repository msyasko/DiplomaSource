using System;
using UnityEngine;

public class UnitySingleton<T> : MonoBehaviour
    where T : Component
{
    private static readonly string instanceName = "UnitySingleton";

    private static T _instance;
    public static T Instance
    {
        get {
            if (_instance == null){
                _instance = FindObjectOfType<T>();
                if (_instance == null){
                    GameObject go = GameObject.Find(instanceName);
                    if (go == null){
                        go = new GameObject();
                        go.name = instanceName;
                    }
                    _instance = go.GetComponent<T>();
                    if (_instance == null) {
                        _instance = go.AddComponent<T>();
                    }
                }
            }
            return _instance;
        }
    }

    public virtual void Awake()
    {
        if (_instance == null) {
            _instance = this as T;
        }
        else { Destroy(this); }
    }
}
