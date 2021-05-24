using UnityEngine;

public class MouseOrbit : MonoBehaviour
{
    public Transform target;
    public float distance = 10.0f;
    public float xSpeed = 250.0f;
    public float ySpeed = 120.0f;
    public float yMinLimit = -20f;
    public float yMaxLimit = 80f;

    public float maxDistance = 100;
    public float minDistance = 0.5f;
    public float dollySpeed = 10;

    private float x;
    private float y;

    private Vector3 startPosP;
    private Quaternion startRotP;

    void OnEnable()
    {
        CalcAngles();
    }

    public void LateUpdate()
    {
        if (!enabled)
            return;

        float dolly = Input.GetAxis("Mouse ScrollWheel");
        float newDist = distance - dolly * dollySpeed * (Mathf.InverseLerp(minDistance*4, maxDistance/4, distance)+0.5f);
        distance = Mathf.Clamp(newDist, minDistance, maxDistance);


        if (target)
        {
            if (Input.GetMouseButton(1))
            {
                x += Input.GetAxis("Mouse X") * xSpeed;
                y -= Input.GetAxis("Mouse Y") * ySpeed;
                y = ClampAngle(y, yMinLimit, yMaxLimit);
            }

            Quaternion rotation = Quaternion.Euler(y, x, 0);
            Vector3 position = (rotation * Vector3.forward * -distance) + target.position;

            transform.rotation = rotation;
            transform.position = position;
        }
    }

    float ClampAngle(float angle, float min, float max)
    {
        if (angle < -360)
            angle += 360;
        if (angle > 360)
            angle -= 360;
        return Mathf.Clamp(angle, min, max);
    }


    [ContextMenu("ResetTransf")]
    private void ResetTransformations()
    {
        transform.position = startPosP;
        transform.rotation = startRotP;

        CalcAngles();
    }

    void Start()
    {
        startPosP = transform.position;
        startRotP = transform.rotation;
    }

    void CalcAngles()
    {
        y = -Vector3.SignedAngle(transform.up, Vector3.up, transform.right);
        Vector3 fwrd = transform.forward;
        fwrd.y = 0;
        x = Vector3.SignedAngle(fwrd, Vector3.forward, Vector3.down);
    }
}