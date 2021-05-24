using System;
using UnityEngine;

public interface IHistoryState : IDisposable
{
    bool Execute();
    bool Undo();
}
