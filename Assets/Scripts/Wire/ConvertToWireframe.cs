using System.Collections.Generic;
using UnityEngine;

public class ConvertToWireframe : MonoBehaviour
{
    void Awake()
    {
        MeshFilter filter = GetComponent<MeshFilter>();

        Mesh origenalMesh = filter.sharedMesh;
        Vector3[] positions = origenalMesh.vertices;
        int[] triangles = origenalMesh.triangles;

        List<int> linesIdexes = new List<int>();

        for (int i = 0; i < triangles.Length-3; i+=3)
        {
            linesIdexes.Add(triangles[i]);
            linesIdexes.Add(triangles[i+1]);

            linesIdexes.Add(triangles[i+1]);
            linesIdexes.Add(triangles[i+2]);

            linesIdexes.Add(triangles[i+2]);
            linesIdexes.Add(triangles[i]);
        }

        Mesh wireMesh = new Mesh();
        wireMesh.vertices = positions;
        wireMesh.SetIndices(linesIdexes.ToArray(), MeshTopology.Lines, 0, false);
        wireMesh.bounds = origenalMesh.bounds;

        filter.mesh = wireMesh;
    }
}