using System;
using UnityEngine;

public class ObjectStateHistory : MonoBehaviour
{
    protected CycleUndoRedoStorage<IHistoryState> _actionHistory;
    [SerializeField] protected int _storageLenght = 5;

    public int Count
    {
        get {
            if (_actionHistory != null) {
                return _actionHistory.Count;
            }
            return 0;
        }
    }

    public int CurrentStateIndex
    {
        get
        {
            if (_actionHistory != null) {
                return _actionHistory.CurrentIndex;
            }
            return -1;
        }
    }

    protected virtual void Awake()
    {
        _actionHistory = new CycleUndoRedoStorage<IHistoryState>(_storageLenght);
    }

    public virtual void AddState(IHistoryState state)
    {
        IHistoryState prevState = _actionHistory.AddReusable(state);
        if (prevState != null) {
            prevState.Dispose();
        }
    }

    public IHistoryState GetState(int index)
    {
        if (_actionHistory == null) { return null; }
        if (index >= Count || index < 0) { return null; }

        return _actionHistory[index];
    }

    public IHistoryState GetCurrentState()
    {
        if (_actionHistory == null) { return null; }
        int curIndex = CurrentStateIndex;
        if (curIndex < 0) { return null; }

        return GetState(curIndex);
    }

    public virtual void Undo()
    {
        _actionHistory.Undo();
    }

    public virtual void Redo()
    {
        _actionHistory.Redo();
    }

    public virtual void Clear()
    {
        _actionHistory.Clear();
    }
}
