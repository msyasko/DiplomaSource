using System;
using System.Collections.Generic;
using ManipulationInstruments;
using UnityEngine;

public class HistoryState_Transform : IHistoryState
{
    public struct Arg
    {
        public Vector3 position;
        public Quaternion rotation;
        public Vector3 scale;

        public Arg (Transform trans)
        {
            if (trans != null)
            {
                position = trans.position;
                rotation = trans.rotation;
                scale = trans.lossyScale;
            }
            else
            {
                position = Vector3.zero;
                rotation = Quaternion.identity;
                scale = Vector3.one;
            }
        }
    }

    private GameObject _go;
    private Arg _from;
    private Arg _to;

    public HistoryState_Transform(GameObject go, Arg from, Arg to)
    {
        _go = go;
        _from = from;
        _to = to;
    }

    public bool Undo()
    {
        if (_go == null) { return false; }

        _go.transform.position = _from.position;
        _go.transform.rotation = _from.rotation;
        _go.transform.localScale = _from.scale;

        return true;
    }

    public bool Execute()
    {
        if (_go == null) { return false; }

        _go.transform.position = _to.position;
        _go.transform.rotation = _to.rotation;
        _go.transform.localScale = _to.scale;

        return true;
    }

    public void Dispose()
    {
        return;
    }
}



public class HistoryState_TransformMultiselect : IHistoryState
{
    public readonly List<TransformSnapshot> oldSnap = new List<TransformSnapshot>();
    public readonly List<TransformSnapshot> newSnap = new List<TransformSnapshot>();

    public void Dispose()
    {
        return;
    }

    public bool Execute()
    {
        foreach (var snap in newSnap)
        {
            snap.transform.position = snap.position;
            snap.transform.rotation= snap.rotation;
            snap.transform.localScale = snap.scale;
        }

        ManipulatorsController.UpdateSelection();

        return true;
    }

    public bool Undo()
    {
        foreach (var snap in oldSnap)
        {
            snap.transform.position = snap.position;
            snap.transform.rotation = snap.rotation;
            snap.transform.localScale = snap.scale;
        }

        ManipulatorsController.UpdateSelection();
        
        return true;
    }

    public struct TransformSnapshot
    {
        public Transform transform;
        public Vector3 position;
        public Quaternion rotation;
        public Vector3 scale;
    }
}

