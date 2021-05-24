using System;

public class CycleUndoRedoStorage<TValue>
    where TValue : IHistoryState
{
    private int _zeroIndex;
    private int _count;
    public int Count {
        get { return _count; }
    }

    private TValue[] _items;
    public TValue[] Items{
        get { return _items; }
    }

    private int _currentIndex = -1;
    public int CurrentIndex
    {
        get { return _currentIndex; }
    }

    public TValue this[int index]
    {
        get
        {
            if (index > _count) { return default (TValue); }
            int ind = (_zeroIndex + index) % _items.Length;
            return _items[ind];
        }
    }

    public CycleUndoRedoStorage(int itemsCount)
    {
        _items = new TValue[itemsCount];
    }

    public void Add(TValue item)
    {
        if (_currentIndex != _count - 1){
            _count = _currentIndex + 1;
        }

        if (_count < _items.Length){
            _count++;
            _currentIndex++;
        }
        else { _zeroIndex++; }

        int index = (_currentIndex + _zeroIndex) % _items.Length;
        _items[index] = item;
    }

    public TValue AddReusable(TValue item)
    {
        if (_currentIndex != _count - 1) {
            _count = _currentIndex + 1;
        }

        if (_count < _items.Length){
            _count++;
            _currentIndex++;
        }
        else { _zeroIndex++; }

        int index = (_currentIndex + _zeroIndex) % _items.Length;
        TValue result = _items[index];
        _items[index] = item;
        return result;
    }

    public TValue Remove(int index)
    {
        TValue result = default(TValue);
        if (_count == 0 || index < 0 || index >= _count) {
            return result;
        }

        result = _items[(_zeroIndex + index) % _items.Length];
        for (int i = index; i < _count - 1; i++)
        {
            int ind = (_zeroIndex + i) % _items.Length;
            int ind2 = (_zeroIndex + i + 1) % _items.Length;
            _items[ind] = _items[ind2];
        }
        _count--;
        if (_currentIndex == _count){
            _currentIndex--;
        }

        return result;
    }

    public bool Undo()
    {
        STEP:
        if (_currentIndex >= 0)
        {
            bool b = _items[(_zeroIndex + _currentIndex) % _items.Length].Undo();
            _currentIndex--;
            if (!b) { goto STEP; }

            return true;
        }
        return false;
    }

    public bool Redo()
    {
        STEP:
        if (_currentIndex < _count - 1)
        {
            _currentIndex++;
            bool b = _items[(_zeroIndex + _currentIndex) % _items.Length].Execute();
            if (!b) { goto STEP; }

            return true;
        }
        return false;
    }

    public void Clear()
    {
        _count = 0;
        _currentIndex = -1;
        _zeroIndex = 0;
        Array.Clear(_items, 0, _count);
    }
}
