using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace Assets.Scrips
{
    public enum Weight { LIGHT, MEDIUM, HEAVY }

    public class TreeNode
    {
        public Vector3 position { get; set; }
        public List<TreeNode> children { get; set; }
        public TreeNode parent { get; set; }
        public Weight weight { get; set; }

        public TreeNode(Vector3 position, TreeNode parent) {
            this.position = position;
            this.parent = parent;
            children = new List<TreeNode>();
        }
    }
}