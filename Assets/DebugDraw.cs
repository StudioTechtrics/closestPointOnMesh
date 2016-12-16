using UnityEngine;
using System.Collections;

public static class DebugDraw {

    public static void DrawBounds(Bounds b, Color c, float t, Transform tr = null)
    {
        Vector3 center = b.center;

        float x = b.extents.x;
        float y = b.extents.y;
        float z = b.extents.z;

        Vector3 ruf = center + new Vector3( x,  y,  z);
        Vector3 rub = center + new Vector3( x,  y, -z);
        Vector3 luf = center + new Vector3(-x,  y,  z);
        Vector3 lub = center + new Vector3(-x,  y, -z);

        Vector3 rdf = center + new Vector3( x, -y,  z);
        Vector3 rdb = center + new Vector3( x, -y, -z);
        Vector3 lfd = center + new Vector3(-x, -y,  z);
        Vector3 lbd = center + new Vector3(-x, -y, -z);

        if (tr != null)
        {
            ruf = tr.TransformPoint(ruf);
            rub = tr.TransformPoint(rub);
            luf = tr.TransformPoint(luf);
            lub = tr.TransformPoint(lub);

            rdf = tr.TransformPoint(rdf);
            rdb = tr.TransformPoint(rdb);
            lfd = tr.TransformPoint(lfd);
            lbd = tr.TransformPoint(lbd);
        }

        Debug.DrawLine(ruf, luf, c, t);
        Debug.DrawLine(ruf, rub, c, t);
        Debug.DrawLine(luf, lub, c, t);
        Debug.DrawLine(rub, lub, c, t);

        Debug.DrawLine(ruf, rdf, c, t);
        Debug.DrawLine(rub, rdb, c, t);
        Debug.DrawLine(luf, lfd, c, t);
        Debug.DrawLine(lub, lbd, c, t);

        Debug.DrawLine(rdf, lfd, c, t);
        Debug.DrawLine(rdf, rdb, c, t);
        Debug.DrawLine(lfd, lbd, c, t);
        Debug.DrawLine(lbd, rdb, c, t);

    }

    public static void DrawMarker(Vector3 position, float size, Color color, float duration, bool depthTest = true)
    {
        Vector3 line1PosA = position + Vector3.up * size * 0.5f;
        Vector3 line1PosB = position - Vector3.up * size * 0.5f;

        Vector3 line2PosA = position + Vector3.right * size * 0.5f;
        Vector3 line2PosB = position - Vector3.right * size * 0.5f;

        Vector3 line3PosA = position + Vector3.forward * size * 0.5f;
        Vector3 line3PosB = position - Vector3.forward * size * 0.5f;

        Debug.DrawLine(line1PosA, line1PosB, color, duration, depthTest);
        Debug.DrawLine(line2PosA, line2PosB, color, duration, depthTest);
        Debug.DrawLine(line3PosA, line3PosB, color, duration, depthTest);
    }

    // Courtesy of robertbu
    public static void DrawPlane(Vector3 position, Vector3 normal, float size, Color color, float duration, bool depthTest = true) 
    {
        Vector3 v3;
 
        if (normal.normalized != Vector3.forward)
            v3 = Vector3.Cross(normal, Vector3.forward).normalized * normal.magnitude;
        else
            v3 = Vector3.Cross(normal, Vector3.up).normalized * normal.magnitude;;
 
        Vector3 corner0 = position + v3 * size;
        Vector3 corner2 = position - v3 * size;
 
        Quaternion q = Quaternion.AngleAxis(90.0f, normal);
        v3 = q * v3;
        Vector3 corner1 = position + v3 * size;
        Vector3 corner3 = position - v3 * size;

        Debug.DrawLine(corner0, corner2, color, duration, depthTest);
        Debug.DrawLine(corner1, corner3, color, duration, depthTest);
        Debug.DrawLine(corner0, corner1, color, duration, depthTest);
        Debug.DrawLine(corner1, corner2, color, duration, depthTest);
        Debug.DrawLine(corner2, corner3, color, duration, depthTest);
        Debug.DrawLine(corner3, corner0, color, duration, depthTest);
        Debug.DrawRay(position, normal * size, color, duration, depthTest);
    }

    public static void DrawVector(Vector3 position, Vector3 direction, float raySize, float markerSize, Color color, float duration, bool depthTest = true)
    {
        Debug.DrawRay(position, direction * raySize, color, 0, false);
        DebugDraw.DrawMarker(position + direction * raySize, markerSize, color, 0, false);
    }

    public static void DrawTriangle(Vector3 a, Vector3 b, Vector3 c, Color color)
    {
        
        Debug.DrawLine(a, b, color);
        Debug.DrawLine(b, c, color);
        Debug.DrawLine(c, a, color);
    }

    public static void DrawTriangle(Vector3 a, Vector3 b, Vector3 c, Color color, Transform t, float duration = 0.1f)
    {
        a = t.TransformPoint(a);
        b = t.TransformPoint(b);
        c = t.TransformPoint(c);

        Debug.DrawLine(a, b, color, duration);
        Debug.DrawLine(b, c, color, duration);
        Debug.DrawLine(c, a, color, duration);
    }

    public static void DrawMesh(Mesh mesh, Color color, Transform t)
    {
        for (int i = 0; i < mesh.triangles.Length; i += 3)
        {
            DrawTriangle(mesh.vertices[mesh.triangles[i]], mesh.vertices[mesh.triangles[i + 1]], mesh.vertices[mesh.triangles[i + 2]], color, t);
        }
    }

    public static Color RandomColor()
    {
        return new Color(Random.value, Random.value, Random.value);
    }
}
