using ManipulationInstruments.Handles;
using UnityEngine;
using UnityEngine.EventSystems;

namespace ManipulationInstruments.UI
{
    public class ToggleButton : MonoBehaviour, IPointerClickHandler
    {
        public GameObject targetObject;


        public void OnPointerClick(PointerEventData eventData)
        {
            targetObject.SetActive(!targetObject.activeSelf);
        }
    }
}