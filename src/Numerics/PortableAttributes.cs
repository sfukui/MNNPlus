﻿#if PORTABLE

using System;

namespace MathNet.Numerics
{
    [AttributeUsage(AttributeTargets.Class | AttributeTargets.Struct)]
    public class SerializableAttribute : Attribute
    {
    }

    [AttributeUsage(AttributeTargets.Constructor | AttributeTargets.Method, AllowMultiple = false, Inherited = false)]
    public class TargetedPatchingOptOutAttribute : Attribute
    {
        public string Reason { get; private set; }
        public TargetedPatchingOptOutAttribute(string reason)
        {
            Reason = reason;
        }
    }

    [AttributeUsage(AttributeTargets.Method, AllowMultiple = false, Inherited = false)]
    public class SpecialNameAttribute : Attribute
    {
    }
}

#endif
