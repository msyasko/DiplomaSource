using UnityEngine;

namespace ManipulationInstruments.Handles
{
    public class AliginUpAxisAlongWorld : MonoBehaviour
    {
        void LateUpdate()
        {
            transform.up = Vector3.up;
        }
    }
}