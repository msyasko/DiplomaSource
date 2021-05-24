using UnityEngine;

public class UIbyStateController : MonoBehaviour
{
    public GameObject calculateSplash;
    public GameObject[] trasformationControls;

    void Awake()
    {
        StateController.ESetState += OnSetState;
    }

    void Destroy()
    {
        StateController.ESetState -= OnSetState;
    }

    private void OnSetState(SceneState state)
    {
        switch (state)
        {
            case SceneState.Setup:
                SetTranformationActivnes(true);
                break;
            case SceneState.Calculation:
                calculateSplash.SetActive(true);
                SetTranformationActivnes(false);
                break;
            case SceneState.Spectate:
                calculateSplash.SetActive(false);
                break;
        }

        void SetTranformationActivnes(bool a)
        {
            foreach (var tc in trasformationControls)
                tc.SetActive(a);
        }
    }
}