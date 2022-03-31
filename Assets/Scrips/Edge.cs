using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace Assets.Scrips
{
    public class Edge
    {
        public Vector3 A { get; set; }
        public Vector3 B { get; set; }
        public float weight { get; set; }

        public Edge(Vector3 A, Vector3 B) {
            this.A = A;
            this.B = B;
            weight = Vector3.Dot(A - B, Vector3.right);
        }
    }
}