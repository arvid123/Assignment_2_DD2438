using UnityEditor;
using UnityEngine;

namespace Assets.Scrips
{
    public class Qnode
    {
        public int i { get; set; }
        public int j { get; set; }
        public Qnode parent { get; set; }

        public Qnode(int i, int j, Qnode parent) {
            this.i = i;
            this.j = j;
            this.parent = parent;
        }
    }
}