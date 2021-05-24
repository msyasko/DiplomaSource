using System.Collections;
using System.Collections.Generic;
using ManipulationInstruments.Handles;
using UnityEngine.EventSystems;
using UnityEngine;

namespace ManipulationInstruments.UI
{
    public class ScaleUpDown : MonoBehaviour, IPointerClickHandler
    {
        //public HandleType handleType;
        public GameObject deactevatedOnClick;

        public string scaleSide;
        private float scaleStep = 0.5f;
        public GameObject cameraObj;

        public void OnPointerClick(PointerEventData eventData)
        {
            if (scaleSide == "Up")
            {
                Vector3 newPos = new Vector3(cameraObj.transform.position.x, cameraObj.transform.position.y - scaleStep, cameraObj.transform.position.z);
                cameraObj.transform.position = newPos;
            }
            else 
            {
                Vector3 newPos = new Vector3(cameraObj.transform.position.x, cameraObj.transform.position.y + scaleStep, cameraObj.transform.position.z);
                cameraObj.transform.position = newPos;
            }
        }
    }
}