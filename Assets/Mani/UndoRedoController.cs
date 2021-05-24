using System;
using UnityEngine;

public class UndoRedoController : MonoBehaviour {

    [SerializeField] private SceneStateHistory _history;
    public SceneStateHistory History
    {
        get { return _history; }
    }

    void Start()
    {
        if (_history == null){ _history = UnitySingleton<SceneStateHistory>.Instance; }
    }

    public void Undo()
    {
        _history.Undo();
    }

    public void Redo()
    {
        _history.Redo();
    }
}
