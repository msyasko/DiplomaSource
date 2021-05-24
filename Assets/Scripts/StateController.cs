using System;
using UnityEngine;

public class StateController : MonoBehaviour
{
    public static event Action<SceneState> ESetState;

    private static SceneState stateS;

    public static SceneState State
    {
        get { return stateS; }
        private set
        {
            stateS = value;
            ESetState?.Invoke(stateS);
        }
    }

    public static void SetState(SceneState state)
    {
        State = state;
    }

    public void SwitchState()
    {
        if (State==SceneState.Setup)
        {
            State = SceneState.Calculation;
        }
        else if(State==SceneState.Spectate)
        {
            State = SceneState.Setup;
        }
    }
}

public enum SceneState
{
    Setup,
    Calculation,
    Spectate
}