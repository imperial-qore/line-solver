import numpy as np
import torch
import torch.nn as nn


class QueueNetworkLearningRNNLayer(nn.Module):
    """RNN layer for learning queueing network parameters.

    Adapted from:
    Garbi, G et al. (2020). Learning Queueing Networks by Recurrent Neural Networks

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """

    def __init__(self, M, R, concurrency):
        super().__init__()
        self.M = M
        self.R = R
        self.concurrency = np.asarray(concurrency, dtype=np.float32)
        self.I = np.eye(M, dtype=np.float32)

        # Learnable parameters
        self.mu = nn.Parameter(self._init_uniform_nonneg((M, 1)))
        P_raw = self._init_uniform_nonneg((M, M - 1))
        P_raw = P_raw / P_raw.sum(dim=1, keepdim=True)
        # Expand to M x M with zeros on diagonal
        P_full = torch.zeros(M, M)
        for i in range(M):
            col = 0
            for j in range(M):
                if i != j:
                    P_full[i, j] = P_raw[i, col]
                    col += 1
        self.P = nn.Parameter(P_full)

        self.hidden_state = None
        self.reset_state()

    def _init_uniform_nonneg(self, sz):
        a = 0.01
        b = 10.0
        return torch.FloatTensor(*sz).uniform_(a, b)

    def reset_state(self):
        self.hidden_state = torch.zeros(self.M + 1)

    def forward(self, X):
        """Forward pass through the RNN layer.

        Args:
            X: tensor of shape (num_time_steps, M, R+1)
               where X[t, :, 0] contains timestamps and X[t, :, 1:] contains queue lengths

        Returns:
            Z: tensor of shape (num_time_steps, M, 2) with [timestamp, predicted_ql]
        """
        num_time_steps = X.shape[0]
        self.reset_state()
        Z = torch.zeros(num_time_steps, self.M, 2)

        concurrency_t = torch.tensor(self.concurrency, dtype=torch.float32)
        I_t = torch.tensor(self.I, dtype=torch.float32)

        for t in range(num_time_steps):
            current_t = X[t, 0, 0]
            old_t = self.hidden_state[0]
            delta_t = current_t - old_t

            pm = torch.abs(self.mu).T * (torch.abs(self.P) - I_t)
            state_ql = self.hidden_state[1:]
            pred = state_ql + (delta_t * torch.min(state_ql, concurrency_t)) @ pm

            if delta_t == 0:
                pred = X[t, :, 1]

            xh_pred = torch.cat([current_t.unsqueeze(0), pred])
            self.hidden_state = xh_pred.detach()

            Z[t, :, 0] = current_t
            Z[t, :, 1] = xh_pred[1:]

        return Z
