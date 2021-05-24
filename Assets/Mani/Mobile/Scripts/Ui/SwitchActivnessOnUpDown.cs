using UnityEngine;
using UnityEngine.EventSystems;

namespace ManipulationInstruments.UI
{
    public class SwitchActivnessOnUpDown:MonoBehaviour, IPointerDownHandler, IPointerUpHandler
    {
        public Behaviour activeOnDown;
        public Behaviour activeOnUp;

        public void OnPointerDown(PointerEventData eventData)
        {
            activeOnDown.enabled = true;
            activeOnUp.enabled = false;
        }

        public void OnPointerUp(PointerEventData eventData)
        {
            activeOnDown.enabled = false;
            activeOnUp.enabled = true;
        }
    }
}