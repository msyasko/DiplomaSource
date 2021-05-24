using UnityEngine;
using UnityEngine.EventSystems;

namespace ManipulationInstruments.Handles
{
    public class ToggleWprldLocal : MonoBehaviour, IPointerClickHandler
    {
        public HandleBase handle;


        public void OnPointerClick(PointerEventData eventData)
        {
            handle.SwithSpace();
        }
    }
}