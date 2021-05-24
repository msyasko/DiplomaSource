using System.Collections.Generic;
using UnityEngine;

public class WireBox : MonoBehaviour
{
    [ContextMenu("Convert")]
    void Awake()
    {
        MeshFilter filter = GetComponent<MeshFilter>();
        List<Vector3> positions = new List<Vector3>();
        List<int> linesIdexes = new List<int>();

        positions.Add(new Vector3(-0.5f, 0.5f, 0.5f));
        positions.Add(new Vector3(0.5f, 0.5f, 0.5f));
        positions.Add(new Vector3(0.5f, 0.5f, -0.5f));
        positions.Add(new Vector3(-0.5f, 0.5f, -0.5f));

        positions.Add(new Vector3(-0.5f, -0.5f, 0.5f));
        positions.Add(new Vector3(0.5f, -0.5f, 0.5f));
        positions.Add(new Vector3(0.5f, -0.5f, -0.5f));
        positions.Add(new Vector3(-0.5f, -0.5f, -0.5f));


        linesIdexes.Add(0);
        linesIdexes.Add(1);
        linesIdexes.Add(1);
        linesIdexes.Add(2);
        linesIdexes.Add(2);
        linesIdexes.Add(3);
        linesIdexes.Add(3);
        linesIdexes.Add(0);

        linesIdexes.Add(4);
        linesIdexes.Add(5);
        linesIdexes.Add(5);
        linesIdexes.Add(6);
        linesIdexes.Add(6);
        linesIdexes.Add(7);
        linesIdexes.Add(7);
        linesIdexes.Add(4);

        linesIdexes.Add(0);
        linesIdexes.Add(4);
        linesIdexes.Add(1);
        linesIdexes.Add(5);
        linesIdexes.Add(2);
        linesIdexes.Add(6);
        linesIdexes.Add(3);
        linesIdexes.Add(7);
        
        Mesh wireMesh = new Mesh();
        wireMesh.vertices = positions.ToArray();
        wireMesh.SetIndices(linesIdexes.ToArray(), MeshTopology.Lines, 0, true);

        filter.mesh = wireMesh;
    } 
}